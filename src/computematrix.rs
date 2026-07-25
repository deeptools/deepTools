use crate::calc::{max_float, mean_float, median_float, min_float, sum_float};
use crate::covcalc::{Bin, Gtfparse, Region, Scalingregions};
use crate::filehandler::{
    bwintervals, chrombounds_from_bw, header_matrix, is_bed_or_gtf, read_bedfile, read_gtffile,
    write_matrix,
};
use itertools::Itertools;
use pyo3::prelude::*;
use pyo3::types::PyList;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
use std::path::Path;

#[pyfunction]
pub fn r_computematrix(
    py: Python,
    mode: &str,                   // reference-point or scale-regions
    regionlis: Py<PyList>,        // python list of region files (bed or gtf)
    bwlis: Py<PyList>,            // python list of bigwig files
    sampleslabel: Py<PyList>,     // python list of sample labels, if empty, use bigwig file names.
    upstream: u32,                // upstream region to consider
    downstream: u32,              // downstream region to consider
    unscaled5prime: u32, // unscaled region 5' of the anchorpoint, only used in scale-regions mode.
    unscaled3prime: u32, // unscaled region 3' of the anchorpoint, only used in scale-regions mode.
    regionbodylength: u32, // length of the region body (after scaling), only used in scale-regions mode.
    binsize: u32,          // binsize to use for the matrix
    missingdatazero: bool, // Encode missing data as 0. Default is false (and will be encoded as NA).
    metagene: bool,        // If set, 'exons' are stitched together to form a metagene
    txnid: &str,           // transcript id to use when parsing GTF file
    exonid: &str,          // exon id to use when parsing GTF file
    txniddesignator: &str, // designator to use when parsing GTF file
    scale: f32,            // scaling factor for writing out values. default is 1.0 (no scaling)
    nanafterend: bool,     // end regions will treated as nans. Default is false.
    skipzeros: bool,       // skip regions with all zeros. Default is false.
    minthresh: f32,        // minimum threshold to keep a region. If not set it will equal 0.0
    maxthresh: f32,        // maximum threshold to keep a region. if not set it will equal 0.0
    averagetypebins: &str, // operation to summarize values over bins. Default is mean.
    sortregions: &str, // either ascend, descend or keep. Default is keep (and ignores sortusing).
    sortusing: &str, // metric to sort on. Either mean median max min sum region_length. Default is mean.
    sortusingsamples: Py<PyList>, // list of samples to sort on. If empty, use all samples.
    referencepoint: &str, // reference point to use. Either TSS, TES or center. Default is TSS. Only used in reference-point mode.
    nproc: usize,         // number of threads.
    verbose: bool,        // verbose output.
    ofile: &str,          // npz file to write to.
    outfilenamematrix: Option<String>, // optional raw matrix tab output file
    outfilesortedregions: Option<String>, // optional sorted/filtered BED output file
    startlabel: Option<String>,  // label for start of region (default "TSS")
    endlabel: Option<String>,    // label for end of region (default "TES")
) -> PyResult<()> {
    // Extract the bed and bigwig files from pyList to Vec.
    let region_files: Vec<String> = regionlis
        .extract(py)
        .expect("Failed to retrieve bed files.");

    let bw_files: Vec<String> = bwlis
        .extract(py)
        .expect("Failed to retrieve bigwig filess.");

    let mut samples_label: Vec<String> = sampleslabel
        .extract(py)
        .expect("Failed to retrieve samples label.");

    let sort_using_samples: Vec<u32> = sortusingsamples
        .extract(py)
        .expect("Failed to retrieve the samples to sort on.");
    // Assert that samples_label equals bw_files, if samples_label is not empty.
    if !samples_label.is_empty() {
        assert_eq!(
            samples_label.len(),
            bw_files.len(),
            "Number of samplelabels do not match number of bigwig files."
        );
    }
    // Assert that sort_using_samples is smaller or equal to bw_files.
    assert!(
        sort_using_samples.len() <= bw_files.len(),
        "Number of samples to sort on is larger than number of bigwig files provided."
    );
    // Get chromosome boundaries from first bigwig file.
    let chromsizes = chrombounds_from_bw(&bw_files.get(0).unwrap());
    // compute number of columns
    let bpsum = &upstream + &downstream + &unscaled5prime + &unscaled3prime + &regionbodylength;

    // Binsize divisibility validation (Feature 11)
    if regionbodylength % binsize != 0 {
        eprintln!("The --regionBodyLength has to be a multiple of --binSize.\nCurrently the values are {} and {} for regionsBodyLength and binSize respectively.", regionbodylength, binsize);
        std::process::exit(1);
    }
    if downstream % binsize != 0 {
        eprintln!("Length of region after the body has to be a multiple of --binSize.\nCurrent value is {}", downstream);
        std::process::exit(1);
    }
    if upstream % binsize != 0 {
        eprintln!("Length of region before the body has to be a multiple of --binSize.\nCurrent value is {}", upstream);
        std::process::exit(1);
    }
    if unscaled5prime % binsize != 0 {
        eprintln!("Length of the unscaled 5 prime region has to be a multiple of --binSize.\nCurrent value is {}", unscaled5prime);
        std::process::exit(1);
    }
    if unscaled3prime % binsize != 0 {
        eprintln!("Length of the unscaled 3 prime region has to be a multiple of --binSize.\nCurrent value is {}", unscaled3prime);
        std::process::exit(1);
    }
    if regionbodylength == 0 && (unscaled5prime > 0 || unscaled3prime > 0) {
        eprintln!("Unscaled 5- and 3-prime regions only make sense with the scale-regions subcommand.");
        std::process::exit(1);
    }

    // Reference-point validation (Feature 12)
    let valid_referencepoint = if mode == "reference-point" {
        if !["TSS", "TES", "center"].contains(&referencepoint) {
            eprintln!("referencepoint must be one of 'TSS', 'TES', or 'center'. Got '{}'", referencepoint);
            std::process::exit(1);
        }
        referencepoint.to_string()
    } else {
        String::new()
    };

    // nanAfterEnd validation (Feature 13): only valid in reference-point mode
    if mode != "reference-point" && nanafterend {
        eprintln!("--nanAfterEnd is only valid in reference-point mode.");
        std::process::exit(1);
    }

    // Get the 'basepaths' of the bed files to use as labels later on (Feature 8: smartLabels).
    let mut regionlabels: Vec<String> = Vec::new();
    for bed in region_files.iter() {
        let entryname = Path::new(bed)
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .into_owned();
        regionlabels.push(entryname);
    }
    if samples_label.is_empty() {
        // no samples labels provided via CLI, retrieve them from bigwig names (Feature 8: smartLabels).
        for bw in bw_files.iter() {
            let entryname = Path::new(bw)
                .file_stem()
                .unwrap()
                .to_string_lossy()
                .into_owned();
            samples_label.push(entryname);
        }
    }
    // Define the scaling regions in a struct
    let scale_regions = Scalingregions {
        upstream: upstream,
        downstream: downstream,
        unscaled5prime: unscaled5prime,
        unscaled3prime: unscaled3prime,
        regionbodylength: regionbodylength,
        binsize: binsize,
        cols_expected: ((bw_files.len() * bpsum as usize) / binsize as usize),
        bpsum: bpsum,
        missingdata_as_zero: missingdatazero,
        scale: scale,
        nan_after_end: nanafterend,
        skipzero: skipzeros,
        minthresh: minthresh,
        maxthresh: maxthresh,
        referencepoint: valid_referencepoint,
        mode: mode.to_string(),
        bwfiles: bw_files.len(),
        avgtype: averagetypebins.to_string(),
        verbose: verbose,
        proc_number: nproc,
        regionlabels: regionlabels,
        bwlabels: samples_label,
        startlabel: startlabel.unwrap_or_else(|| "TSS".to_string()),
        endlabel: endlabel.unwrap_or_else(|| "TES".to_string()),
    };
    let gtfparse = Gtfparse {
        metagene: metagene,
        txnid: txnid.to_string(),
        exonid: exonid.to_string(),
        txniddesignator: txniddesignator.to_string(),
    };
    if verbose {
        println!("Region files: {:?}", &region_files);
        println!("Bigwig files: {:?}", &bw_files);
        println!("Samples labels: {:?}", scale_regions.bwlabels);
        println!("Sort using samples: {:?}", &sort_using_samples);
    }
    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();

    // Parse regions from bed files. Note that we retain the name of the bed file (in case there are more then 1)
    // Additionaly, score and strand are also retained, if it's a 3-column bed file we just fill in '.'
    let mut regions: Vec<Region> = Vec::new();
    let mut regionsizes: HashMap<String, u32> = HashMap::new();
    region_files.iter()
        .map(|r| {
            let ftype = is_bed_or_gtf(r);

            match ftype.as_str() {
                "gtf" => read_gtffile(r, &gtfparse, chromsizes.keys().collect()),
                "bed" => read_bedfile(r, metagene, chromsizes.keys().collect()),
                _ => panic!("Only .bed and .gtf files are allowed (as determined by the number of columns). File = {}", ftype),
            }
        })
        .for_each(|(reg, regsize)| {
            regions.extend(reg);
            regionsizes.insert(regsize.0, regsize.1);
        });
    // Define slop regions, which contain the actual 'bins' to query the bigwig files.
    let slopregions = pool.install(|| {
        regions
            .par_iter()
            .map(|region| slop_region(&region, &scale_regions, &chromsizes))
            .collect::<Vec<_>>()
    });

    // Discriminate between reference-point and scale-regions mode.

    let matrix: Vec<Vec<f32>> = pool.install(|| {
        bw_files
            .par_iter()
            .map(|i| bwintervals(&i, &regions, &slopregions, &scale_regions))
            .reduce(
                || vec![vec![]; regions.len()],
                |mut acc, vec_of_vecs| {
                    for (i, inner_vec) in vec_of_vecs.into_iter().enumerate() {
                        acc[i].extend(inner_vec);
                    }
                    acc
                },
            )
    });
    matrix_dump(
        sortregions,
        sortusing,
        sort_using_samples,
        regions,
        matrix,
        scale_regions,
        regionsizes,
        ofile,
        outfilenamematrix,
        outfilesortedregions,
        verbose,
    );

    Ok(())
}

fn slop_region(
    region: &Region,
    scale_regions: &Scalingregions,
    chromsizes: &HashMap<String, u32>,
) -> Vec<Bin> {
    // Idea is to create a vector Bins (Conbin or Catbin) which encodes start and end of every bin (binsize passed by computeMatrix).
    // Catbin takes care of the situation where one needs metagenes, and thus multiple start/end per bin are possible.
    // The number of columns is predetermined
    // Note that the before / after could mean that we run out of chromosome.
    // Invalid regions (later to be encoded as NA or 0), will be pushed as (0,0) tuples.
    // Note that if nan_after_end is set to true, we will push (0,0) tuples after the end of the region.

    // Get the chromosome end for a specific region, and assert that the region stays within the chromosome boundary.
    // Note that only a right check is needed, as positions are u32.
    // Note that we know &region.chrom is inside chromsizes already, since this filtering is done at the region reading stage.
    let chromend: u32 = *chromsizes.get(&region.chrom).unwrap();
    region.assert_end(chromend);
    region.get_anchor_bins(scale_regions, chromend)
}

#[allow(unused_mut)]
#[allow(unused_assignments)]
fn matrix_dump(
    sortregions: &str,
    sortusing: &str,
    sort_using_samples: Vec<u32>,
    regions: Vec<Region>,
    matrix: Vec<Vec<f32>>,
    scale_regions: Scalingregions,
    regionsizes: HashMap<String, u32>,
    ofile: &str,
    outfilenamematrix: Option<String>,
    outfilesortedregions: Option<String>,
    verbose: bool,
) {
    // Takes a pre-computed matrix, resorts it if requested, and writes it to file.
    // Resort the matrix, if this is requested.
    if sortregions != "keep" {
        if verbose {
            println!(
                "Sorting output matrix with settings: sortRegions: {}, sortUsing {}",
                sortregions, sortusing
            );
        }
        // If sortusingsamples is set, we need a vector to subset the columns of interest
        let mut cols_of_interest: Vec<usize> = Vec::new();
        if !sort_using_samples.is_empty() {
            let cols_per_sample = scale_regions.cols_expected / scale_regions.bwfiles;
            for sample_ix in sort_using_samples.iter() {
                // Note that sort_using_samples is assumed to be 1-index. Hence we need to subtract 1.
                let start = (sample_ix - 1) * cols_per_sample as u32;
                let end = start + cols_per_sample as u32;
                cols_of_interest.extend(start as usize..end as usize);
            }
        }
        let mut regionslices: Vec<(usize, usize)> = Vec::new();
        let mut rstart: usize = 0;
        let mut rend: usize = 0;
        let mut lastregion = &regions[0].name;
        for (ix, region) in regions.iter().enumerate() {
            if region.name != *lastregion {
                regionslices.push((rstart, rend));
                rstart = ix;
                rend = ix;
                lastregion = &region.name;
            }
            rend = ix;
        }
        regionslices.push((rstart, rend));
        let mut sortedix: Vec<usize>;
        if sortusing == "region_length" {
            if !sort_using_samples.is_empty() && verbose {
                println!("Sort using samples is set ({:?}), but is not used when sorting on region_length. It is thus ignored.", sort_using_samples);
            }
            sortedix = regionslices
                .iter()
                .flat_map(|(start, end)| {
                    let rslice = &regions[*start..*end+1];
                    let tix = rslice
                        .iter()
                        .enumerate()
                        .map(|(ix, region)| {
                            (ix + *start, region.regionlength)
                        })
                        .collect::<Vec<_>>()
                        .iter()
                        .sorted_by(|ix, metric| ix.1.partial_cmp(&metric.1).unwrap())
                        .map(|(ix, _)| *ix)
                        .collect::<Vec<usize>>();
                    match sortregions {
                        "ascend" => tix,
                        "descend" => tix.into_iter().rev().collect(),
                        _ => panic!("If sortRegions is not keep, it should be either ascend or descend. Not {}", sortregions),
                    }
                })
                .collect();
        } else {
            sortedix = regionslices
                .iter()
                .flat_map(|(start, end)| {
                    let rslice = &matrix[*start..*end+1];
                    let tix = rslice
                        .iter()
                        .enumerate()
                        .map(|(ix, vals)| {
                            let subset: Vec<_> = if cols_of_interest.is_empty() {
                                vals.iter().collect()
                            } else {
                                cols_of_interest
                                    .iter()
                                    .filter_map(|&index| vals.get(index))  // `vec.get(index)` returns Option<&T>
                                    .collect()
                            };
                            let metric = match sortusing {
                                "mean" => mean_float(subset),
                                "median" => median_float(subset),
                                "max" => max_float(subset),
                                "min" => min_float(subset),
                                "sum" => sum_float(subset),
                                _ => panic!("Sortusing should be either mean, median, max, min, sum or region_length. Not {}", sortusing),
                            };
                            (ix + *start, metric)
                        })
                        .collect::<Vec<_>>()
                        .iter()
                        .sorted_by(|ix, metric| ix.1.partial_cmp(&metric.1).unwrap())
                        .map(|(ix, _)| *ix)
                        .collect::<Vec<usize>>();
                    match sortregions {
                        "ascend" => tix,
                        "descend" => tix.into_iter().rev().collect(),
                        _ => panic!("If sortRegions is not keep, it should be either ascend or descend. Not {}", sortregions),
                    }
                })
                .collect();
        }
        // assert sorted ix length == matrix length == regions length
        assert_eq!(
            sortedix.len(),
            matrix.len(),
            "Length of sorted indices does not match matrix length: {} != {}",
            sortedix.len(),
            matrix.len()
        );
        assert_eq!(
            sortedix.len(),
            regions.len(),
            "Length of sorted indices does not match regions length: {} ! = {}",
            sortedix.len(),
            regions.len()
        );

        // Reorder matrix & regions
        let sortedmatrix: Vec<Vec<f32>> = sortedix.iter().map(|ix| matrix[*ix].clone()).collect();
        let sortedregions: Vec<Region> = sortedix.into_iter().map(|ix| regions[ix].clone()).collect();
        write_matrix(
            header_matrix(&scale_regions, &regionsizes, sortregions, sortusing),
            sortedmatrix.clone(),
            ofile,
            sortedregions.clone(),
            &scale_regions,
        );

        // Feature 6: outFileNameMatrix - write raw matrix tab file
        if let Some(ref matrix_file) = outfilenamematrix {
            crate::filehandler::write_matrix_values(matrix_file, &sortedmatrix, &scale_regions, &regionsizes);
        }

        // Feature 7: outFileSortedRegions - write sorted BED file
        if let Some(ref bed_file) = outfilesortedregions {
            crate::filehandler::write_sorted_regions_bed(bed_file, &sortedregions, &scale_regions, &regionsizes);
        }
    } else {
        write_matrix(
            header_matrix(&scale_regions, &regionsizes, sortregions, sortusing),
            matrix.clone(),
            ofile,
            regions.clone(),
            &scale_regions,
        );

        // Feature 6: outFileNameMatrix - write raw matrix tab file
        if let Some(ref matrix_file) = outfilenamematrix {
            crate::filehandler::write_matrix_values(matrix_file, &matrix, &scale_regions, &regionsizes);
        }

        // Feature 7: outFileSortedRegions - write regions BED file
        if let Some(ref bed_file) = outfilesortedregions {
            crate::filehandler::write_sorted_regions_bed(bed_file, &regions, &scale_regions, &regionsizes);
        }
    }
}
