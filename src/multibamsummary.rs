use crate::calc::deseq_scalefactors;
use crate::covcalc::{bam_pileup, parse_regions, region_divider, TempZip};
use crate::covcalc::{Gtfparse, Region};
use crate::filehandler::{
    bam_ispaired, chrombounds_from_bam, is_bed_or_gtf, read_bedfile, read_gtffile,
};
use crate::filtering::Alignmentfilters;
use itertools::multiunzip;
use ndarray::Array2;
use ndarray_npy::NpzWriter;
use pyo3::prelude::*;
use pyo3::types::PyList;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};

#[pyfunction]
pub fn r_mbams(
    py: Python,
    // required parameters
    mode: &str, // either bins or BED-file
    bam_files: Py<PyList>,
    ofile: &str,
    // additional output
    out_raw_counts: &str,
    scaling_factors: &str,
    // optional parameters
    labels: Py<PyList>,
    mut binsize: u32,
    distance_between_bins: u32,
    nproc: usize,
    bed_file: Py<PyList>,
    supregion: &str,
    blacklist: &str,
    verbose: bool,
    extendreads: bool,
    extendreadslen: u32,
    centerreads: bool,
    samflaginclude: u16,   // sam flag include
    samflagexclude: u16,   // sam flag exclude
    minfraglen: u32,       // minimum fragment length.
    maxfraglen: u32,       // maximum fragment length.
    minmappingquality: u8, // minimum mapping quality.
    metagene: bool,        // metagene mode or not.
    txnid: &str,           // transcript id to use when parsing GTF file.
    exonid: &str,          // exon id to use when parsing GTF file.
    txniddesignator: &str, // designator to use when parsing GTF file.
) -> PyResult<()> {
    let bamfiles: Vec<String> = bam_files
        .extract(py)
        .expect("Failed to retrieve bam files.");

    let bamlabels: Vec<String> = labels.extract(py).expect("Failed to retrieve labels.");

    let bedfiles: Vec<String> = bed_file.extract(py).expect("Failed to retrieve bedfiles.");

    let ignorechr: Vec<String> = Vec::new();

    let max_len = bamlabels.iter().map(|s| s.len()).max().unwrap_or(0);
    let bamlabels_arr: Array2<u8> = Array2::from_shape_fn((bamlabels.len(), max_len), |(i, j)| {
        bamlabels[i].as_bytes().get(j).copied().unwrap_or(0)
    });

    // Get paired-end information
    let ispe = bamfiles.iter().map(|x| bam_ispaired(x)).collect::<Vec<_>>();

    // zip through ispe and bamfiles
    if verbose {
        for (_ispe, _bf) in ispe.iter().zip(bamfiles.iter()) {
            println!("Sample: {} is-paired: {}", _bf, _ispe);
        }
    }

    let mut regions: Vec<Region> = Vec::new();
    let mut gene_mode = false;
    let mut blacklistregions: Option<Vec<Region>> = None;
    if mode == "BED-file" {
        if verbose {
            println!("BED file mode. with files: {:?}", bedfiles);
        }
        if supregion != "None" {
            println!("Region supplied in BED-file mode. The region will be ignored.");
        }
        let gtfparse = Gtfparse {
            metagene: metagene,
            txnid: txnid.to_string(),
            exonid: exonid.to_string(),
            txniddesignator: txniddesignator.to_string(),
        };
        // From the first bamfile, get the chromosome sizes. Only the chromosome names are needed, to make sure no invalid regions are included later on.
        let chromsizes = chrombounds_from_bam(bamfiles.iter().map(|x| x.as_str()).collect());

        binsize = 1;
        // Parse regions from bed files. Note that we retain the name of the bed file (in case there are more then 1)
        // Additionaly, score and strand are also retained, if it's a 3-column bed file we just fill in '.'
        let mut regionsizes: HashMap<String, u32> = HashMap::new();
        bedfiles.iter()
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
        gene_mode = true;

        // If there is a blacklist, read it.
        if blacklist != "None" {
            // Check if it's a bed or gtf file
            let isbed = is_bed_or_gtf(blacklist);
            match isbed.as_str() {
                "gtf" => panic!("Error: Please provide a bed file for the blacklist."),
                "bed" => {
                    let (bls, _) =
                        read_bedfile(&blacklist.to_string(), false, chromsizes.keys().collect());
                    blacklistregions = Some(bls);
                }
                _ => panic!("Error: Cannot determine filetype of blacklist file."),
            }
        }
    } else {
        if verbose {
            println!(
                "BINS mode. with binsize: {}, distance between bins: {}",
                binsize, distance_between_bins
            );
        }
        let (parsedregions, chromsizes) =
            parse_regions(supregion, bamfiles.iter().map(|x| x.as_str()).collect());
        regions = parsedregions;
        // If there is a blacklist, read it.
        if blacklist != "None" {
            // Check if it's a bed or gtf file
            let isbed = is_bed_or_gtf(blacklist);
            match isbed.as_str() {
                "gtf" => panic!("Error: Please provide a bed file for the blacklist."),
                "bed" => {
                    let (bls, _) =
                        read_bedfile(&blacklist.to_string(), false, chromsizes.keys().collect());
                    blacklistregions = Some(bls);
                }
                _ => panic!("Error: Cannot determine filetype of blacklist file."),
            }
        }
    }

    let filters = Alignmentfilters::new(
        blacklistregions,
        Some(minmappingquality),
        Some(samflaginclude),
        Some(samflagexclude),
        Some(minfraglen),
        Some(maxfraglen),
        None, // No MNase mode.
        None, // No offset
        None, // No strand filtering.
        Some(extendreads),
        Some(extendreadslen),
        Some(centerreads),
    );

    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();

    // Zip together bamfiles and ispe into a vec of tuples.
    let bampfiles: Vec<_> = bamfiles.into_iter().zip(ispe.into_iter()).collect();
    let bam_ispe_filter: Vec<(String, bool, Alignmentfilters)> = bampfiles.into_iter()
        .map(|(bamfile, ispe)| {
            let mut filter = filters.clone();
            if filter.extendreads && filter.extendreadslen == 0 && ispe {
                // We need a pass over the bamfile already to get the mean fragment length.
                filter.set_extendreadslen(&bamfile, nproc, &regions);
                if filter.extendreadslen == 0 {
                    panic!("Error: No fragment length found for read extension. Please provide a valid fragment length.");
                }
                if verbose {
                    println!("fragment length for read extension set as: {} for {}", filter.extendreadslen, bamfile.to_string());
                }
            }
            (bamfile, ispe, filter)
        })
        .collect();
    // Divide up the regions into regionBlocks
    let regionblocks = region_divider(&regions);

    assert!(regionblocks.len() > 0, "No regions to process. Exiting.");
    if verbose {
        println!(
            "Regions divided into {} parallel blocks",
            regionblocks.len()
        );
        println!("Start coverage calculation");
    }

    let covcalcs: Vec<_> = pool.install(|| {
        bam_ispe_filter
            .par_iter()
            .map(|(bamfile, ispe, filter)| {
                let (bg, _mapped, _unmapped, _readlen, _fraglen) = regionblocks
                    .par_iter()
                    .map(|i| {
                        bam_pileup(
                            bamfile, &i, &binsize, &ispe, &ignorechr, filter, false, gene_mode,
                            false, 0,
                        )
                    })
                    .reduce(
                        || (vec![], 0, 0, vec![], vec![]),
                        |(mut _bg, mut _mapped, mut _unmapped, mut _readlen, mut _fraglen),
                         (bg, mapped, unmapped, readlen, fraglen)| {
                            _bg.extend(bg);
                            _readlen.extend(readlen);
                            _fraglen.extend(fraglen);
                            _mapped += mapped;
                            _unmapped += unmapped;
                            (_bg, _mapped, _unmapped, _readlen, _fraglen)
                        },
                    );
                bg
            })
            .collect()
    });
    if verbose {
        println!("Coverage calculation done");
        println!("Define output file");
    }

    // Collate the coverage files into a matrix.
    let its: Vec<_> = covcalcs.iter().map(|x| x.into_iter()).collect();
    let zips = TempZip { iterators: its };
    if verbose {
        println!("Start iterating through temp coverage files and create output npy.");
    }
    let zips_vec: Vec<_> = zips.collect();
    println!(" Length of ziperators = {}", zips_vec.len());

    let matvec: Vec<_> = pool.install(|| {
        let _m: Vec<_> = zips_vec
            .par_iter()
            .flat_map(|c| {
                let readers: Vec<_> = c
                    .par_iter()
                    .map(|x| BufReader::new(File::open(x).unwrap()).lines())
                    .collect();
                let mut _matvec: Vec<Vec<f32>> = Vec::new();
                let mut _regions: Vec<(String, String, String)> = Vec::new();
                for mut _l in (TempZip { iterators: readers }) {
                    // unwrap all lines in _l
                    let lines: Vec<_> = _l
                        .par_iter_mut()
                        .map(|x| x.as_mut().unwrap())
                        .map(|x| x.split('\t').collect())
                        .map(|x: Vec<&str>| {
                            (
                                x[0].to_string(),
                                x[1].to_string(),
                                x[2].to_string(),
                                x[3].parse::<f32>().unwrap(),
                            )
                        })
                        .collect();
                    let counts = lines
                        .par_iter()
                        .map(|x| (x.3 * 100.0).round() / 100.0)
                        .collect::<Vec<_>>();
                    let regions: (String, String, String) =
                        (lines[0].0.clone(), lines[0].1.clone(), lines[0].2.clone());
                    _matvec.push(counts);
                    _regions.push(regions);
                }
                (_matvec, _regions)
            })
            .collect();
        _m
    });

    // Seperate the matrices and regions
    let (matvec, regions): (Vec<_>, Vec<_>) = multiunzip(matvec);

    // Write out the count files, if appropriate

    if out_raw_counts != "None" {
        if verbose {
            println!("Writing raw counts to disk.");
        }
        let mut cfile = BufWriter::new(File::create(out_raw_counts).unwrap());
        // Write the header to the file.
        let mut headstr = String::new();
        headstr.push_str("#\'chr\'\t\'start\'\t\'end\'");
        for label in bamlabels.iter() {
            headstr.push_str(&format!("\t\'{}\'", label));
        }
        writeln!(cfile, "{}", headstr).unwrap();
        let outlines: Vec<String> = pool.install(|| {
            regions
                .par_iter()
                .zip(matvec.par_iter())
                .map(|(region, counts)| {
                    let mut outstr = String::new();
                    outstr.push_str(&format!("{}\t{}\t{}", region.0, region.1, region.2));
                    for count in counts.iter() {
                        outstr.push_str(&format!("\t{}", count));
                    }
                    outstr
                })
                .collect()
        });
        for line in outlines {
            writeln!(cfile, "{}", line).unwrap();
        }
    }

    // Create 2darray from matvec
    let matarr: Array2<f32> = Array2::from_shape_vec(
        (matvec.len(), matvec[0].len()),
        matvec.into_iter().flatten().collect(),
    )
    .unwrap();

    // If scalefactors are required, calc and save them now.
    if scaling_factors != "None" {
        if verbose {
            println!("Calculating scale factors.");
        }
        let sf = deseq_scalefactors(&matarr);
        // save scalefactors to file
        let mut sf_file = File::create(scaling_factors).unwrap();
        writeln!(sf_file, "Sample\tscalingFactor").unwrap();
        for (sf, label) in sf.iter().zip(bamlabels.iter()) {
            writeln!(sf_file, "{}\t{}", label, sf).unwrap();
        }
    }

    if verbose {
        println!("Writing matrix to disk.");
    }
    let mut npz = NpzWriter::new_compressed(File::create(ofile).unwrap());
    npz.add_array("matrix", &matarr).unwrap();
    npz.add_array("labels", &bamlabels_arr).unwrap();
    npz.finish().unwrap();
    if verbose {
        println!("Matrix written.");
    }
    Ok(())
}
