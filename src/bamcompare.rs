use pyo3::prelude::*;
use pyo3::types::PyList;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::io::prelude::*;
use std::io::{BufReader};
use std::fs::File;
use itertools::Itertools;
use bigtools::{Value};
use crate::filehandler::{bam_ispaired, write_covfile};
use crate::covcalc::{bam_pileup, parse_regions, Alignmentfilters, region_divider};
use crate::normalization::scale_factor_bamcompare;
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
    operation: &str,
    pseudocount: f32,
    // filtering options
    ignoreduplicates: bool,
    minmappingquality: u8, // 
    samflaginclude: u16,
    samflagexclude: u16,
    minfraglen: u32,
    maxfraglen: u32,
    nproc: usize,
    _ignorechr: Py<PyList>,
    binsize: u32,
    regions: Vec<(String, u32, u32)>,
    verbose: bool
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
    // Set alignment filters
    let filters = Alignmentfilters {
        minmappingquality: minmappingquality,
        samflaginclude: samflaginclude,
        samflagexclude: samflagexclude,
        minfraglen: minfraglen,
        maxfraglen: maxfraglen
    };

    // Parse regions & calculate coverage. Note that 
    let (regions, chromsizes)  = parse_regions(&regions, vec![bamifile1, bamifile2]);
    let regionblocks = region_divider(&regions);

    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();
    
    // Set up the bam files in a Vec.
    let bamfiles = vec![(bamifile1, ispe1), (bamifile2, ispe2)];

    let covcalcs: Vec<ParsedBamFile> = pool.install(|| {
        bamfiles.par_iter()
            .map(|(bamfile, ispe)| {
                let (bg, mapped, unmapped, readlen, fraglen) = regionblocks.par_iter()
                    .map(|i| bam_pileup(bamfile, &i, &binsize, &ispe, &ignorechr, &filters, false))
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

    // Calculate scale factors.
    let sf = scale_factor_bamcompare(scalefactorsmethod, covcalcs[0].mapped, covcalcs[1].mapped, binsize, effective_genome_size, norm);
    println!("scale factor1 = {}, scale factor2 = {}", sf.0, sf.1);
    // Create output stream
    let mut chrom = "".to_string();
    let lines = covcalcs[0].bg.iter().zip(covcalcs[1].bg.iter()).flat_map(
        |(t1, t2)| {
            let reader1 = BufReader::new(File::open(t1).unwrap()).lines();
            let reader2 = BufReader::new(File::open(t2).unwrap()).lines();

            reader1.zip(reader2).map(
                |(l1, l2)| {
                    let l1 = l1.unwrap();
                    let l2 = l2.unwrap();
                    let fields1: Vec<&str> = l1.split('\t').collect();
                    let fields2: Vec<&str> = l2.split('\t').collect();
        
                    let chrom1: String = fields1[0].to_string();
                    let chrom2: String = fields2[0].to_string();
                    let start1: u32 = fields1[1].parse().unwrap();
                    let start2: u32 = fields2[1].parse().unwrap();
                    let end1: u32 = fields1[2].parse().unwrap();
                    let end2: u32 = fields2[2].parse().unwrap();
        
                    // Assert the regions are equal.
                    assert_eq!(chrom1, chrom2);
                    assert_eq!(start1, start2);
                    assert_eq!(end1, end2);
        
                    // Calculate the coverage.
                    let cov1: f32 = fields1[3].parse().unwrap();
                    let cov2: f32 = fields2[3].parse().unwrap();
                    let cov = calc_ratio(cov1, cov2, &sf.0, &sf.1, &pseudocount, operation);
        
                    (chrom1, Value { start: start1, end: end1, value: cov })
                }).coalesce(|p, c| {
                if p.1.value == c.1.value {
                    Ok((p.0, Value {start: p.1.start, end: c.1.end, value: p.1.value}))
                } else {
                    Err((p, c))
                }
            })
        }
    );
    write_covfile(lines, ofile, ofiletype, chromsizes);
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