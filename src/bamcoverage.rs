use pyo3::prelude::*;
use pyo3::types::PyList;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use core::panic;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use bigtools::Value;
use crate::covcalc::{bam_pileup, parse_regions, region_divider, Region};
use crate::filtering::Alignmentfilters;
use crate::filehandler::{bam_ispaired, write_covfile, is_bed_or_gtf, read_bedfile};
use crate::normalization::scale_factor;
use crate::calc::median;


#[pyfunction]
pub fn r_bamcoverage(
    // input and output
    bamifile: &str, // input bamfile
    ofile: &str, // output file
    ofiletype: &str, // output file type, bedgraph or bigwig
    // norm options
    norm: &str, // normalization mode RPKM, CPM, BPM, RPGC
    effectivegenomesize: u64, // default is 0, when not set. 
    scalefactor: f32, // default 1.0
    // processing options
    mnase: bool,
    _offset: Py<PyList>, // list of 2 [offset 5', offset 3'], if no offset is required we have [0, 0]
    extendreads: bool, // true for extension of reads
    extendreadslen: u32, // if extendreads is set, and SE, this length is used for extension.
    centerreads: bool, // to center the reads or not.
    filterrnastrand: &str, // forward, reverse or 'None'
    blacklist: &str, // path to blacklist filename, or 'None'
    _ignorechr: Py<PyList>, // list of chromosomes to ignore. Is empty if none.
    skipnoncovregions: bool,
    _smoothlength: u32, // 0 = no smoothing, else it's a strictly larger then binsize
    binsize: u32,
    // filtering options
    minmappingquality: u8, // 
    samflaginclude: u16,
    samflagexclude: u16,
    mut minfraglen: u32, // default 0 -> no filter
    mut maxfraglen: u32, // defualt 0 -> no filter
    nproc: usize,
    supregion: &str,
    verbose: bool,
    collapse: bool,
) -> PyResult<()> {
    let mut offset: (i32, i32) = (0, 0);
    let mut ignorechr: Vec<String> = Vec::new();
    Python::with_gil(|py| {
        let offsetv: Vec<i32> = _offset.extract(py).expect("Failed to retrieve offset.");
        if offsetv.len() == 2 {
            offset = (offsetv[0], offsetv[1]);
        } else {
            panic!("Error: Offset should be a list of 2. Received: {:?}", offsetv);
        }
        ignorechr = _ignorechr.extract(py).expect("Failed to retrieve ignorechr.");
    });
    let ispe = bam_ispaired(bamifile);

    // If mnase, library should be PE !
    if mnase && !ispe {
        panic!("Error: MNase-seq requires paired-end data.");
    }
    if norm == "RPGC" && effectivegenomesize == 0 {
        panic!("Error: Effective genome size is required for RPGC normalization.");
    }
    if norm != "None" && scalefactor != 1.0 {
        println!("Warning: You have set a normalization option ({}), but also a scale factor. Only the scale factor will be used", norm);
    }
    if mnase && offset != (0, 0) {
        println!("Warning: Both MNase and offset are set. The offset will be ignored !");
    }
    if offset != (0, 0) {
        if offset.1 > 0 && offset.1 < offset.0 {
            panic!("Right side offset cannot be smaller than the left side offset.");
        }
    }
    if extendreads && extendreadslen == 0 && !ispe {
        panic!("Error: Extendreads is set, but not to a specific length (and library is single end). Specify the length with the --extendReads parameter.");
    }
    if centerreads && !extendreads {
        println!("Warning: Centerreads is set, but extendreads is not. Centering will do nothing in this case.");
    }

    if verbose {
        println!("Chromosomes to ignore for normalization: {:?}", ignorechr);
    }

    // if mnase, set the min / max fragment lengths if these are not set.
    if mnase {
        if minfraglen == 0 {
            minfraglen = 130;
        } else if minfraglen != 130{
            println!("Note that MNase mode is set, but minfraglen is set at {}. Recommended is 130.", minfraglen);
        }
        if maxfraglen == 0 {
            maxfraglen = 200;
        } else if maxfraglen != 200 {
            println!("Note that MNase mode is set, but maxfraglen is set at {}. Recommended is 200.", maxfraglen);
        }
        if binsize != 1 {
            println!("Note that MNase mode is set, but binsize is set at {}. Recommended is 1.", binsize);
        }
    }
    if verbose {
        println!("Sample: {} is-paired: {}", bamifile, ispe);
    }
    // Parse regions & calculate coverage
    let (regions, chromsizes)  = parse_regions(supregion, vec![bamifile]);
    let regionblocks = region_divider(&regions);

    // If there is a blacklist, read it.
    let mut backlistregions: Option<Vec<Region>> = None;
    if blacklist != "None" {
        // Check if it's a bed or gtf file
        let isbed = is_bed_or_gtf(blacklist);
        match isbed.as_str() {
            "gtf" => panic!("Error: Please provide a bed file for the blacklist."),
            "bed" => {
                let (bls, _) = read_bedfile(&blacklist.to_string(), false, chromsizes.keys().collect());
                backlistregions = Some(bls);
            },
            _ => panic!("Error: Cannot determine filetype of blacklist file.")
        }
    }
    // Set alignment filters
    let mut filters = Alignmentfilters::new(
        backlistregions,
        Some(minmappingquality),
        Some(samflaginclude),
        Some(samflagexclude),
        Some(minfraglen),
        Some(maxfraglen),
        Some(mnase),
        Some(offset),
        Some(filterrnastrand.to_string()),
        Some(extendreads),
        Some(extendreadslen),
        Some(centerreads),
    );
    // If extendreads, extendreadslen & ispe, we need a pass over the bamfile already now to get the mean fraglen.
    if filters.extendreads && filters.extendreadslen == 0 && ispe {
        // We need a pass over the bamfile already to get the mean fragment length.
        filters.set_extendreadslen(bamifile, nproc, &regions);
        if verbose {
            println!("fragment length for read extension set as: {}", filters.extendreadslen);
        }
    }

    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();
    let (bg, mapped, _unmapped, readlen, fraglen) = pool.install(|| {
        regionblocks.par_iter()
            .map(|i| bam_pileup(bamifile, &i, &binsize, &ispe, &ignorechr, &filters, collapse, false, true))
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
            )
    });

    let mut readlen = median(readlen);
    if filters.extendreads {
        if verbose {
            println!("extend reads option on, overriding readlen from {} to {}", readlen, filters.extendreadslen);
        }
        readlen = filters.extendreadslen as f32;
    }
    let fraglen = median(fraglen);
    if verbose {
        println!("Read stats with ignorechr: {:?}", ignorechr);
        println!("Mapped: {} Unmapped: {}", mapped, _unmapped);
        println!("Readlen: {}, Fraglen: {}", readlen, fraglen);
    }

    let sf = scale_factor(
        norm, 
        mapped,
        binsize,
        effectivegenomesize,
        readlen,
        fraglen,
        scalefactor,
        &verbose
    );

    // Create output stream
    let lines = bg.into_iter().flat_map(
        |bg| {
            let reader = BufReader::new(File::open(bg).unwrap());
            reader.lines().filter_map(
                |l| {
                    let l = l.unwrap();
                    let fields: Vec<&str> = l.split('\t').collect();
                    if skipnoncovregions && fields[3] == "0" {
                        None
                    } else {
                        Some(
                            (fields[0].to_string(), Value {
                                start: fields[1].parse::<u32>().unwrap(),
                                end: fields[2].parse::<u32>().unwrap(),
                                value: (fields[3].parse::<f32>().unwrap() * sf * 100.0).round() / 100.0,
                            })
                        )
                    }
                }
            )
        }
    );
    if verbose {
        println!("Writing output to: {}", ofile);
    }

    write_covfile(lines, ofile, ofiletype, chromsizes);
    Ok(())
}
