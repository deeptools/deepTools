use pyo3::prelude::*;
use pyo3::types::PyList;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::io::prelude::*;
use std::io::{BufReader};
use std::fs::File;
use itertools::Itertools;
use bigtools::{Value};
use crate::filehandler::{bam_ispaired, write_covfile, is_bed_or_gtf, read_bedfile};
use crate::covcalc::{bam_pileup, parse_regions, TempZip, region_divider, Region};
use crate::filtering::Alignmentfilters;
use crate::normalization::{scale_factor_bamcompare};
use crate::calc::{median, calc_ratio};
use tempfile::{TempPath};

#[pyfunction]
pub fn r_bamcompare(
    // input and output
    bamifile1: &str, // input bamfile 1
    bamifile2: &str, // input bamfile 2
    ofile: &str, // output file
    ofiletype: &str, // ouput file type, bedgraph or bigwig
    // norm options
    norm: &str,
    effective_genome_size: u64,
    scalefactorsmethod: &str,
    given_sf1: f32,
    given_sf2: f32,
    operation: &str,
    pseudocount1: f32,
    pseudocount2: f32,
    // filtering options
    extendreads: bool, // if 0, no extension
    extendreadslen: u32, // length of extension (0 if PE or if not extending)
    centerreads: bool,
    blacklist: &str, // path to blacklist filename, or 'None'
    minmappingquality: u8, // 
    samflaginclude: u16,
    samflagexclude: u16,
    minfraglen: u32,
    maxfraglen: u32,
    nproc: usize,
    _ignorechr: Py<PyList>,
    skip_non_covered_regions: bool,
    skip_zero_over_zero: bool,
    binsize: u32,
    supregion: &str,
    verbose: bool,
    collapse: bool
) -> PyResult<()> {
    let ispe1 = bam_ispaired(bamifile1);
    let ispe2 = bam_ispaired(bamifile2);

    if verbose {
        println!("Sample1: {} is-paired: {}", bamifile1, ispe1);
        println!("Sample2: {} is-paired: {}", bamifile2, ispe2);
    }
    let mut ignorechr: Vec<String> = Vec::new();
    Python::with_gil(|py| {
        ignorechr = _ignorechr.extract(py).expect("Failed to retrieve ignorechr.");
    });


    // Parse regions & calculate coverage. Note that 
    let (regions, chromsizes)  = parse_regions(supregion, vec![bamifile1, bamifile2]);
    let regionblocks = region_divider(&regions);

    // If there is a blacklist, read it.
    let mut blacklistregions: Option<Vec<Region>> = None;
    if blacklist != "None" {
        // Check if it's a bed or gtf file
        let isbed = is_bed_or_gtf(blacklist);
        match isbed.as_str() {
            "gtf" => panic!("Error: Please provide a bed file for the blacklist."),
            "bed" => {
                let (bls, _) = read_bedfile(&blacklist.to_string(), false, chromsizes.keys().collect());
                blacklistregions = Some(bls);
            },
            _ => panic!("Error: Cannot determine filetype of blacklist file.")
        }
    }
    // 
    // Set alignment filters
    let mut filter1 = Alignmentfilters::new(
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
    let mut filter2 = filter1.clone();
    if filter1.extendreads && filter1.extendreadslen == 0 && !ispe1 {
        panic!("Error: Extendreads is set, but not to a specific length (and library is single end). Specify the length with the --extendReads parameter.");
    }
    if filter2.extendreads && filter2.extendreadslen == 0 && !ispe2 {
        panic!("Error: Extendreads is set, but not to a specific length (and library is single end). Specify the length with the --extendReads parameter.");
    }
    if filter1.extendreads && filter1.extendreadslen == 0 && ispe1 {
        // We need a pass over the bamfile already to get the mean fragment length.
        filter1.set_extendreadslen(bamifile1, nproc, &regions);
        if verbose {
            println!("fragment length for read extension set as: {} for bamfile 1", filter1.extendreadslen);
        }
    }
    if filter2.extendreads && filter2.extendreadslen == 0 && ispe2 {
        // We need a pass over the bamfile already to get the mean fragment length.
        filter2.set_extendreadslen(bamifile2, nproc, &regions);
        if verbose {
            println!("fragment length for read extension set as: {} for bamfile 2", filter2.extendreadslen);
        }
    }

    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();
    
    // Set up the bam files in a Vec.
    let bamfiles: Vec<(&str, bool, &Alignmentfilters)> = vec![(bamifile1, ispe1, &filter1), (bamifile2, ispe2, &filter2)];

    let mut covcalcs: Vec<ParsedBamFile> = pool.install(|| {
        bamfiles.par_iter()
            .map(|(bamfile, ispe, alfilter)| {
                let (bg, mapped, unmapped, readlen, fraglen) = regionblocks.par_iter()
                    .map(|i| bam_pileup(bamfile, &i, &binsize, &ispe, &ignorechr, alfilter , false, false, true))
                    .reduce(
                        || (vec![], 0, 0, vec![], vec![]),
                        |(mut _bg, mut _mapped, mut _unmapped, mut _readlen, mut _fraglen), (bg, mapped, unmapped, readlen, fraglen)| {
                            _bg.extend(bg);
                            _readlen.extend(readlen);
                            _fraglen.extend(fraglen);
                            _mapped += mapped;
                            _unmapped += unmapped;
                            (_bg, _mapped, _unmapped, _readlen, _fraglen)
                        }
                    );
                ParsedBamFile {
                    bamfile: bamfile,
                    ispe: *ispe,
                    bg: bg,
                    mapped: mapped,
                    unmapped: unmapped,
                    readlen: median(readlen),
                    fraglen: median(fraglen)
                }
            })
        .collect()
    });
    // Print out some stats if verbose
    if verbose {
        println!("bamfile\tPE\tmapped\tunmapped\tmed_readlen\tmed_fraglen");
        println!("{}\t{}\t{}\t{}\t{}\t{}", covcalcs[0].bamfile, covcalcs[0].ispe, covcalcs[0].mapped, covcalcs[0].unmapped, covcalcs[0].readlen, covcalcs[0].fraglen);
        println!("{}\t{}\t{}\t{}\t{}\t{}", covcalcs[1].bamfile, covcalcs[1].ispe, covcalcs[1].mapped, covcalcs[1].unmapped, covcalcs[1].readlen, covcalcs[1].fraglen);
    }
    // Calculate scale factors.
    
    let mut sf = scale_factor_bamcompare(
        scalefactorsmethod,
        covcalcs[0].mapped, covcalcs[1].mapped,
        binsize, effective_genome_size, norm,
        covcalcs[0].readlen, covcalcs[1].readlen,
        covcalcs[0].fraglen, covcalcs[1].fraglen,
    );
    if given_sf1 != 0.0 || given_sf2 != 0.0 {
        if verbose {
            println!("Using given scale factors: {} and {}", given_sf1, given_sf2);
        }
        sf = (given_sf1, given_sf2);
    } else if verbose {
        println!("scale factor1 = {}, scale factor2 = {}", sf.0, sf.1);
    }
    
    // Extract both vecs of TempPaths into a single vector
    let its = vec![
        covcalcs[0].bg.drain(..).collect::<Vec<_>>(),
        covcalcs[1].bg.drain(..).collect::<Vec<_>>()
    ];
    let its: Vec<_> = its.iter().map(|x| x.into_iter()).collect();
    let zips = TempZip { iterators: its };
    let zips_vec: Vec<_> = zips.collect();

    if collapse {
        let lines = zips_vec
            .into_iter()
            .flat_map(|c| {
                let readers: Vec<_> = c.into_iter().map(|x| BufReader::new(File::open(x).unwrap()).lines()).collect();
                let temp_zip = TempZip { iterators: readers };
                temp_zip.into_iter().filter_map(|mut _l| {
                    let lines: Vec<_> = _l
                        .iter_mut()
                        .map(|x| x.as_mut().unwrap())
                        .map(|x| x.split('\t').collect())
                        .map(|x: Vec<&str>| (x[0].to_string(), x[1].parse::<u32>().unwrap(), x[2].parse::<u32>().unwrap(), x[3].parse::<f32>().unwrap()))
                        .collect();
                    assert_eq!(lines.len(), 2);
                    assert_eq!(lines[0].0, lines[1].0, "Error: Chromosome mismatch in bam files. {} != {}", lines[0].0, lines[1].0);
                    assert_eq!(lines[0].1, lines[1].1, "Error: Start position mismatch in bam files. {} != {}", lines[0].1, lines[1].1);
                    assert_eq!(lines[0].2, lines[1].2, "Error: End position mismatch in bam files. {} != {}", lines[0].2, lines[1].2);
                    // Calculate the coverage.
                    if skip_zero_over_zero && lines[0].3 == 0.0 && lines[1].3 == 0.0 {
                        return None;
                    } else if skip_non_covered_regions && lines[0].3 == 0.0 {
                        return None;
                    } else if skip_non_covered_regions && lines[1].3 == 0.0 {
                        return None;
                    } else {
                        let cov = calc_ratio(lines[0].3, lines[1].3, &sf.0, &sf.1, &pseudocount1, &pseudocount2, operation);
                        Some((lines[0].0.clone(), Value { start: lines[0].1, end: lines[0].2, value: cov }))
                    }
                }).coalesce(|p, c| {
                    if p.1.value == c.1.value && p.0 == c.0 {
                        Ok((p.0, Value {start: p.1.start, end: c.1.end, value: p.1.value}))
                    } else {
                        Err((p, c))
                    }
                })
        });
        write_covfile(lines, ofile, ofiletype, chromsizes);
    } else {
        let lines = zips_vec
            .into_iter()
            .flat_map(|c| {
                let readers: Vec<_> = c.into_iter().map(|x| BufReader::new(File::open(x).unwrap()).lines()).collect();
                let temp_zip = TempZip { iterators: readers };
                temp_zip.into_iter().filter_map(|mut _l| {
                    let lines: Vec<_> = _l
                        .iter_mut()
                        .map(|x| x.as_mut().unwrap())
                        .map(|x| x.split('\t').collect())
                        .map(|x: Vec<&str>| (x[0].to_string(), x[1].parse::<u32>().unwrap(), x[2].parse::<u32>().unwrap(), x[3].parse::<f32>().unwrap()))
                        .collect();
                    assert_eq!(lines.len(), 2);
                    assert_eq!(lines[0].0, lines[1].0, "Error: Chromosome mismatch in bam files. {} != {}", lines[0].0, lines[1].0);
                    assert_eq!(lines[0].1, lines[1].1, "Error: Start position mismatch in bam files. {} != {}", lines[0].1, lines[1].1);
                    assert_eq!(lines[0].2, lines[1].2, "Error: End position mismatch in bam files. {} != {}", lines[0].2, lines[1].2);
                    // Calculate the coverage.
                    if skip_zero_over_zero && lines[0].3 == 0.0 && lines[1].3 == 0.0 {
                        return None;
                    } else if skip_non_covered_regions && lines[0].3 == 0.0 {
                        return None;
                    } else if skip_non_covered_regions && lines[1].3 == 0.0 {
                        return None;
                    } else {
                        let cov = calc_ratio(lines[0].3, lines[1].3, &sf.0, &sf.1, &pseudocount1, &pseudocount2, operation);
                        Some((lines[0].0.clone(), Value { start: lines[0].1, end: lines[0].2, value: cov }))
                    }
                })
        });
        write_covfile(lines, ofile, ofiletype, chromsizes);
    }
    Ok(())
}

pub struct ParsedBamFile<'a> {
    pub bamfile: &'a str,
    pub ispe: bool,
    pub bg: Vec<TempPath>,
    pub mapped: u32,
    pub unmapped: u32,
    pub readlen: f32,
    pub fraglen: f32
}