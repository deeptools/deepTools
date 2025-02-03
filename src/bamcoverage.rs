use pyo3::prelude::*;
use pyo3::types::PyList;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use bigtools::Value;
use crate::covcalc::{bam_pileup, parse_regions, Alignmentfilters, region_divider, Region};
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
    _offset: Py<PyList>, // list of max 2 [offset 5', offset 3'], if no offset is required we have [1, -1]
    _extendreads: u32, // if 0, no extension
    _centerreads: bool,
    _filterrnastrand: &str, // forward, reverse or 'None'
    blacklist: &str, // path to blacklist filename, or 'None'
    _ignorechr: Py<PyList>, // list of chromosomes to ignore. Is empty if none.
    _skipnoncovregions: bool,
    _smoothlength: u32, // 0 = no smoothing, else it's a strictly larger then binsize
    binsize: u32,
    // filtering options
    _ignoreduplicates: bool,
    minmappingquality: u8, // 
    samflaginclude: u16,
    samflagexclude: u16,
    minfraglen: u32,
    maxfraglen: u32,
    nproc: usize,
    supregion: &str,
    verbose: bool
) -> PyResult<()> {
    let mut offset: Vec<i32> = Vec::new();
    let mut ignorechr: Vec<String> = Vec::new();
    Python::with_gil(|py| {
        offset = _offset.extract(py).expect("Failed to retrieve offset.");
        ignorechr = _ignorechr.extract(py).expect("Failed to retrieve ignorechr.");
    });
    let ispe = bam_ispaired(bamifile);
    validate_args(
        &ispe, &mnase,
        &norm, &effectivegenomesize,
        &scalefactor,
        &offset, &ignorechr, &verbose
    );
    // Set alignment filters
    let filters = Alignmentfilters {
        minmappingquality: minmappingquality,
        samflaginclude: samflaginclude,
        samflagexclude: samflagexclude,
        minfraglen: minfraglen,
        maxfraglen: maxfraglen
    };
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

    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();
    let (bg, mapped, _unmapped, readlen, fraglen) = pool.install(|| {
        regionblocks.par_iter()
            .map(|i| bam_pileup(bamifile, &i, &binsize, &ispe, &ignorechr, &filters, true, false, &backlistregions))
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
    let readlen = median(readlen);
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
            reader.lines().map(
                |l| {
                    let l = l.unwrap();
                    let fields: Vec<&str> = l.split('\t').collect();
                    let chrom: String = fields[0].to_string();
                    let start: u32 = fields[1].parse().unwrap();
                    let end: u32 = fields[2].parse().unwrap();
                    let cov: f32 = fields[3].parse().unwrap();
                    (chrom, Value {start: start, end: end, value: cov * sf})
                }
            )
        }
    );

    write_covfile(lines, ofile, ofiletype, chromsizes);
    Ok(())
}

fn validate_args(
    ispe: &bool,
    mnase: &bool,
    norm: &str,
    effectivegenomesize: &u64,
    scalefactor: &f32,
    _offset: &Vec<i32>,
    ignorechr: &Vec<String>,
    verbose: &bool
) {
    // If mnase, library should be PE !
    if *mnase && !ispe {
        panic!("Error: MNase-seq requires paired-end data.");
    }
    if norm == "RPGC" && *effectivegenomesize == 0 {
        panic!("Error: Effective genome size is required for RPGC normalization.");
    }
    if norm != "None" && *scalefactor != 1.0 {
        println!("Warning: You have set a normalization option ({}), but also a scale factor. Only the scale factor will be used", norm);
    }
    if *verbose {
        println!("Chromosomes to ignore for normalization: {:?}", ignorechr);
    }
}