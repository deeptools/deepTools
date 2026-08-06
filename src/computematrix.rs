use crate::calc::{max_float, mean_float, median_float, min_float, sum_float};
use crate::covcalc::{Bin, Gtfparse, Region, Scalingregions};
use crate::filehandler::{
    bwintervals, chrombounds_from_bw, header_matrix, is_bed_or_gtf, read_bedfile, read_gtffile,
    write_matrix,
};
use crate::filtering::BlacklistIndex;
use itertools::Itertools;
use pyo3::prelude::*;
use pyo3::types::PyList;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

#[pyfunction]
pub fn r_computematrix(
    py: Python,
    mode: &str,                           // reference-point or scale-regions
    regionlis: Py<PyList>,                // python list of region files (bed or gtf)
    bwlis: Py<PyList>,                    // python list of bigwig files
    sampleslabel: Py<PyList>, // python list of sample labels, if empty, use bigwig file names.
    upstream: u32,            // upstream region to consider
    downstream: u32,          // downstream region to consider
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
    blacklist: &str,       // path to blacklist BED file, or "None"
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
    startlabel: Option<String>, // label for start of region (default "TSS")
    endlabel: Option<String>, // label for end of region (default "TES")
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

    // region / binsize divisibility validation
    assert!(
        regionbodylength % binsize == 0,
        "The --regionBodyLength has to be a multiple of --binSize."
    );
    assert!(
        downstream % binsize == 0,
        "Length of region after the body has to be a multiple of --binSize."
    );
    assert!(
        upstream % binsize == 0,
        "Length of region before the body has to be a multiple of --binSize."
    );
    assert!(
        unscaled5prime % binsize == 0,
        "Length of the unscaled 5 prime region has to be a multiple of --binSize."
    );
    assert!(
        unscaled3prime % binsize == 0,
        "Length of the unscaled 3 prime region has to be a multiple of --binSize."
    );

    // Reference-point validation
    let valid_referencepoint = if mode == "reference-point" {
        assert!(
            ["TSS", "TES", "center"].contains(&referencepoint),
            "referencepoint must be one of 'TSS', 'TES', or 'center'. Got '{}'",
            referencepoint
        );
        referencepoint.to_string()
    } else {
        String::new()
    };

    // nanAfterEnd validation only valid in reference-point mode
    if mode != "reference-point" {
        assert!(
            !nanafterend,
            "--nanAfterEnd is only valid in reference-point mode."
        )
    }

    // If there is a blacklist, read it and build an index.
    let blacklist_index: Option<Arc<BlacklistIndex>> = if blacklist != "none" {
        let isbed = is_bed_or_gtf(blacklist);
        match isbed.as_str() {
            "gtf" => panic!("Error: Please provide a bed file for the blacklist."),
            "bed" => {
                let (bls, _) =
                    read_bedfile(&blacklist.to_string(), false, chromsizes.keys().collect());
                let idx = BlacklistIndex::from_regions(&bls);
                Some(Arc::new(idx))
            }
            _ => panic!("Error: Cannot determine filetype of blacklist file."),
        }
    } else {
        None
    };

    // Get the 'basepaths' of the bed files to use as labels later on
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
        // no samples labels provided via CLI, retrieve them from bigwig names
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
    let mut regionslices: Vec<(usize, usize)> = Vec::new();
    let mut cursor: usize = 0;

    // Pre-allocate capacity for regions to avoid multiple reallocations
    let total_regions_estimate = region_files
        .iter()
        .map(|r| {
            let ftype = is_bed_or_gtf(r);
            match ftype.as_str() {
                "gtf" => 1000, // rough estimate for GTF files
                "bed" => 1000, // rough estimate for BED files
                _ => 0,
            }
        })
        .sum::<usize>();

    regions.reserve(total_regions_estimate);

    // Collect region data in parallel and then merge
    let region_data: Vec<(Vec<Region>, (String, u32))> = region_files.par_iter()
        .map(|r| {
            let ftype = is_bed_or_gtf(r);

            match ftype.as_str() {
                "gtf" => read_gtffile(r, &gtfparse, chromsizes.keys().collect()),
                "bed" => read_bedfile(r, metagene, chromsizes.keys().collect()),
                _ => panic!("Only .bed and .gtf files are allowed (as determined by the number of columns). File = {}", ftype),
            }
        })
        .collect();

    // Merge regions and region sizes
    // Track which group (region file) each region belongs to, by index.
    let mut group_of_region: Vec<usize> = Vec::with_capacity(cursor);
    for (group_idx, (reg, regsize)) in region_data.into_iter().enumerate() {
        let n = reg.len();
        regions.extend(reg);
        regionsizes.insert(regsize.0.clone(), regsize.1);

        if n > 0 {
            regionslices.push((cursor, cursor + n - 1));
            group_of_region.extend(std::iter::repeat(group_idx).take(n));
            cursor += n;
        }
    }

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
            .map(|i| {
                bwintervals(
                    &i,
                    &regions,
                    &slopregions,
                    &scale_regions,
                    blacklist_index.as_deref(),
                )
            })
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
        group_of_region,
        matrix,
        scale_regions,
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

fn should_skip_row(row: &[f32], scale_regions: &Scalingregions) -> bool {
    // Rows with all nan values are skipped.
    if row.iter().all(|&x| x.is_nan()) {
        return true;
    }

    if scale_regions.skipzero && row.iter().all(|&x| x == 0.0) {
        return true;
    }
    if scale_regions.minthresh != 0.0
        && !row.iter().any(|&x| x.is_nan())
        && row.iter().any(|&x| x <= scale_regions.minthresh)
    {
        return true;
    }
    if scale_regions.maxthresh != 0.0
        && !row.iter().any(|&x| x.is_nan())
        && row.iter().any(|&x| x >= scale_regions.maxthresh)
    {
        return true;
    }
    false
}

fn recompute_regionsizes(
    group_of_region: &[usize],
    scale_regions: &Scalingregions,
) -> (HashMap<String, u32>, Vec<(usize, usize)>) {
    let mut filtered_regionsizes: HashMap<String, u32> = HashMap::new();
    for label in &scale_regions.regionlabels {
        filtered_regionsizes.insert(label.clone(), 0);
    }
    for &group_idx in group_of_region {
        *filtered_regionsizes
            .get_mut(&scale_regions.regionlabels[group_idx])
            .unwrap() += 1;
    }

    // Build filtered regionslices in label order, matching the sorting logic expectations.
    let mut slices: Vec<(usize, usize)> = Vec::new();
    let mut cursor = 0usize;
    for label in &scale_regions.regionlabels {
        let count = *filtered_regionsizes.get(label).unwrap() as usize;
        if count > 0 {
            slices.push((cursor, cursor + count - 1));
            cursor += count;
        }
    }

    (filtered_regionsizes, slices)
}

fn matrix_dump(
    sortregions: &str,
    sortusing: &str,
    sort_using_samples: Vec<u32>,
    regions: Vec<Region>,
    group_of_region: Vec<usize>,
    matrix: Vec<Vec<f32>>,
    scale_regions: Scalingregions,
    ofile: &str,
    outfilenamematrix: Option<String>,
    outfilesortedregions: Option<String>,
    verbose: bool,
) {
    // Filter rows first (before sorting/header generation) so group_boundaries are correct.
    let keep_indices: Vec<usize> = (0..matrix.len())
        .filter(|&i| !should_skip_row(&matrix[i], &scale_regions))
        .collect();

    if verbose {
        println!(
            "Filtering: {} regions -> {} regions",
            regions.len(),
            keep_indices.len()
        );
    }

    let filtered_matrix: Vec<Vec<f32>> = keep_indices.iter().map(|&i| matrix[i].clone()).collect();
    let filtered_regions: Vec<Region> = keep_indices.iter().map(|&i| regions[i].clone()).collect();
    let filtered_group_of_region: Vec<usize> =
        keep_indices.iter().map(|&i| group_of_region[i]).collect();
    let (filtered_regionsizes, filtered_regionslices) =
        recompute_regionsizes(&filtered_group_of_region, &scale_regions);

    if sortregions != "keep" {
        if verbose {
            println!(
                "Sorting output matrix with settings: sortRegions: {}, sortUsing {}",
                sortregions, sortusing
            );
        }
        let mut cols_of_interest: Vec<usize> = Vec::new();
        if !sort_using_samples.is_empty() {
            let cols_per_sample = scale_regions.cols_expected / scale_regions.bwfiles;
            for sample_ix in sort_using_samples.iter() {
                let start = (sample_ix - 1) * cols_per_sample as u32;
                let end = start + cols_per_sample as u32;
                cols_of_interest.extend(start as usize..end as usize);
            }
        }
        if verbose {
            println!(
                "regionslices: {} slices for {} regions",
                filtered_regionslices.len(),
                filtered_regions.len()
            );
        }
        let mut sortedix: Vec<usize>;
        if sortregions == "no" {
            if verbose && (sortusing != "mean" || !sort_using_samples.is_empty()) {
                println!(
                    "sortRegions is 'no': sorting by genomic coordinates (chrom, start, end). sortUsing/sortUsingSamples are ignored."
                );
            }
            sortedix = filtered_regionslices
                .iter()
                .flat_map(|(start, end)| {
                    let rslice = &filtered_regions[*start..*end + 1];
                    rslice
                        .iter()
                        .enumerate()
                        .map(|(ix, region)| (ix + *start, region))
                        .sorted_by(|a, b| {
                            a.1.chrom
                                .cmp(&b.1.chrom)
                                .then_with(|| a.1.start.start_key().cmp(&b.1.start.start_key()))
                                .then_with(|| a.1.end.end_key().cmp(&b.1.end.end_key()))
                        })
                        .map(|(ix, _)| ix)
                        .collect::<Vec<usize>>()
                })
                .collect();
        } else if sortusing == "region_length" {
            if !sort_using_samples.is_empty() && verbose {
                println!(
                    "Sort using samples is set ({:?}), but is not used when sorting on region_length. It is thus ignored.",
                    sort_using_samples
                );
            }
            sortedix = filtered_regionslices
                .iter()
                .flat_map(|(start, end)| {
                    let rslice = &filtered_regions[*start..*end+1];
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
            sortedix = filtered_regionslices
                .iter()
                .flat_map(|(start, end)| {
                    let rslice = &filtered_matrix[*start..*end+1];
                    let tix = rslice
                        .iter()
                        .enumerate()
                        .map(|(ix, vals)| {
                            let subset: Vec<_> = if cols_of_interest.is_empty() {
                                vals.iter().collect()
                            } else {
                                cols_of_interest
                                    .iter()
                                    .filter_map(|&index| vals.get(index))
                                    .collect()
                            };
                            let metric = match sortusing {
                                "mean" => mean_float(&subset),
                                "median" => median_float(&subset),
                                "max" => max_float(&subset),
                                "min" => min_float(&subset),
                                "sum" => sum_float(&subset),
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
        assert_eq!(
            sortedix.len(),
            filtered_matrix.len(),
            "Length of sorted indices does not match matrix length: {} != {}",
            sortedix.len(),
            filtered_matrix.len()
        );
        assert_eq!(
            sortedix.len(),
            filtered_regions.len(),
            "Length of sorted indices does not match regions length: {} ! = {}",
            sortedix.len(),
            filtered_regions.len()
        );

        let sortedmatrix: Vec<Vec<f32>> = sortedix
            .iter()
            .map(|ix| filtered_matrix[*ix].clone())
            .collect();
        let sortedregions: Vec<Region> = sortedix
            .into_iter()
            .map(|ix| filtered_regions[ix].clone())
            .collect();
        write_matrix(
            header_matrix(
                &scale_regions,
                &filtered_regionsizes,
                sortregions,
                sortusing,
            ),
            &sortedmatrix,
            ofile,
            &sortedregions,
            &scale_regions,
        );

        if let Some(ref matrix_file) = outfilenamematrix {
            crate::filehandler::write_matrix_values(
                matrix_file,
                &sortedmatrix,
                &scale_regions,
                &filtered_regionsizes,
            );
        }

        if let Some(ref bed_file) = outfilesortedregions {
            crate::filehandler::write_sorted_regions_bed(
                bed_file,
                &sortedregions,
                &scale_regions,
                &filtered_regionsizes,
            );
        }
    } else {
        write_matrix(
            header_matrix(
                &scale_regions,
                &filtered_regionsizes,
                sortregions,
                sortusing,
            ),
            &filtered_matrix,
            ofile,
            &filtered_regions,
            &scale_regions,
        );

        if let Some(ref matrix_file) = outfilenamematrix {
            crate::filehandler::write_matrix_values(
                matrix_file,
                &filtered_matrix,
                &scale_regions,
                &filtered_regionsizes,
            );
        }

        if let Some(ref bed_file) = outfilesortedregions {
            crate::filehandler::write_sorted_regions_bed(
                bed_file,
                &filtered_regions,
                &scale_regions,
                &filtered_regionsizes,
            );
        }
    }
}
