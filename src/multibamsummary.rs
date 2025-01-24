use pyo3::prelude::*;
use pyo3::types::PyList;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use itertools::{multizip, multiunzip, izip};
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use ndarray::{Array1, Array2, stack, Axis};
use std::io;
use ndarray_npy::NpzWriter;
use std::borrow::Cow;
use std::collections::HashMap;
use std::path::Path; 
use std::sync::{Arc, Mutex};
use crate::covcalc::{bam_pileup, parse_regions, Alignmentfilters, region_divider};
use crate::filehandler::{bam_ispaired, read_bedfile, read_gtffile, chrombounds_from_bam, is_bed_or_gtf};
use crate::calc::{median, calc_ratio, deseq_scalefactors};
use crate::bamcompare::ParsedBamFile;
use crate::normalization::scale_factor_bamcompare;
use crate::covcalc::{Region, Gtfparse};

#[pyfunction]
pub fn r_mbams(
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
    sup_regions: Vec<(String, u32, u32)>,
    _blacklist: &str,
    verbose: bool,
    _extend_reads: u32,
    _center_reads: bool,
    sam_flag_incl: u16, // sam flag include
    sam_flag_excl: u16, // sam flag exclude
    min_fragment_length: u32, // minimum fragment length.
    max_fragment_length: u32, // maximum fragment length.
    min_mapping_quality: u8, // minimum mapping quality.
    metagene: bool, // metagene mode or not.
    txnid: &str, // transcript id to use when parsing GTF file.
    exonid: &str, // exon id to use when parsing GTF file.
    txniddesignator: &str, // designator to use when parsing GTF file.
) -> PyResult<()> {
    let mut bamfiles: Vec<String> = Vec::new();
    let mut bamlabels: Vec<String> = Vec::new();
    let mut bedfiles: Vec<String> = Vec::new();
    let mut ignorechr: Vec<String> = Vec::new();
    Python::with_gil(|py| {
        bamfiles = bam_files.extract(py).expect("Failed to retrieve bam files.");
        bamlabels = labels.extract(py).expect("Failed to retrieve labels.");
        bedfiles = bed_file.extract(py).expect("Failed to retrieve bedfiles.");
    });

    let max_len = bamlabels.iter().map(|s| s.len()).max().unwrap_or(0);
    let bamlabels_arr: Array2<u8> = Array2::from_shape_fn((bamlabels.len(), max_len), |(i, j)| {
        bamlabels[i].as_bytes().get(j).copied().unwrap_or(0)
    });

    // Get paired-end information
    let ispe = bamfiles.iter()
        .map(|x| bam_ispaired(x))
        .collect::<Vec<_>>();

    // zip through ispe and bamfiles
    if verbose {
        for (_ispe, _bf) in ispe.iter().zip(bamfiles.iter()) {
            println!("Sample: {} is-paired: {}", _bf, _ispe);
        }
    }
    
    let filters: Alignmentfilters = Alignmentfilters {
        minmappingquality: min_mapping_quality,
        samflaginclude: sam_flag_incl,
        samflagexclude: sam_flag_excl,
        minfraglen: min_fragment_length,
        maxfraglen: max_fragment_length
    };
    
    let mut regions: Vec<Region> = Vec::new();

    if mode == "BED-file" {
        if verbose {
            println!("BED file mode. with files: {:?}", bedfiles);
        }
        let gtfparse = Gtfparse {
            metagene: metagene,
            txnid: txnid.to_string(),
            exonid: exonid.to_string(),
            txniddesignator: txniddesignator.to_string(),
        };
        // From the first bamfile, get the chromosome sizes. Only the chromosome names are needed, to make sure no invalid regions are included later on.
        let chromsizes = chrombounds_from_bam(bamfiles.get(0).unwrap());

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
    } else {
        if verbose {
            println!("BINS mode. with binsize: {}", binsize);
        }
        let (parsedregions, chromsizes) = parse_regions(&sup_regions, bamfiles.get(0).unwrap());
        regions = parsedregions;
    }

    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();    
    
    // Zip together bamfiles and ispe into a vec of tuples.
    let bampfiles: Vec<_> = bamfiles.into_iter().zip(ispe.into_iter()).collect();

    // Divide up the regions into regionBlocks
    let regionblocks = region_divider(&regions);
    
    assert!(regionblocks.len() > 0, "No regions to process. Exiting.");
    if verbose {
        println!("Regions divided into {} parallel blocks", regionblocks.len());
        println!("Start coverage calculation");
    }

    let covcalcs: Vec<_> = pool.install(|| {
        bampfiles.par_iter()
            .map(|(bamfile, ispe)| {
                let (bg, _mapped, _unmapped, _readlen, _fraglen) = regionblocks.par_iter()
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
    let matvec: Array2<f32> = pool.install(|| {
        let matvecs: Vec<Array1<f32>> = zips_vec
            .par_iter()
            .flat_map(|c| {
                let readers: Vec<_> = c.iter().map(|x| BufReader::new(File::open(x).unwrap()).lines()).collect();
                let mut _matvec: Vec<Array1<f32>> = Vec::new();
                for mut _l in (TempZip { iterators: readers }) {
                    // unwrap all lines in _l
                    let lines: Vec<_> = _l
                        .iter_mut()
                        .map(|x| x.as_mut().unwrap())
                        .map(|x| x.split('\t').collect())
                        .map(|x: Vec<&str> | ( x[0].to_string(), x[1].parse::<u32>().unwrap(), x[2].parse::<u32>().unwrap(), x[3].parse::<f32>().unwrap() ) )
                        .collect();
                    let arrline = Array1::from(lines.iter().map(|x| x.3).collect::<Vec<_>>());
                    _matvec.push(arrline);
                }
                _matvec
            })
            .collect();
        stack(Axis(0), &matvecs.iter().map(|x| x.view()).collect::<Vec<_>>()).unwrap()
    });

    // If scalefactors are required, calc and save them now.
    if scaling_factors != "None" {
        let sf = deseq_scalefactors(&matvec);
        // save scalefactors to file
        let mut sf_file = File::create(scaling_factors).unwrap();
        writeln!(sf_file, "Sample\tscalingFactor").unwrap();
        for (sf, label) in sf.iter().zip(bamlabels.iter()) {
            writeln!(sf_file, "{}\t{}", label, sf).unwrap();
        }
    }
    let mut npz = NpzWriter::new(File::create(ofile).unwrap());
    npz.add_array("matrix", &matvec).unwrap();
    npz.add_array("labels", &bamlabels_arr).unwrap();
    npz.finish().unwrap();
    if verbose {
        println!("Matrix written.");
    }

    // If raw counts are required, write them to file.
    if out_raw_counts != "None" {
        let mut cfile = File::create(out_raw_counts).unwrap();
        // Create the header
        let mut headstr = String::new();
        headstr.push_str("#\'chr\'\t\'start\'\t\'end\'");
        for label in bamlabels.iter() {
            headstr.push_str(&format!("\t\'{}\'", label));
        }
        writeln!(cfile, "{}", headstr).unwrap();
        for (i, region) in regions.iter().enumerate() {
            let mut line = String::new();
            line.push_str(&format!("{}\t{}\t{}", region.chrom, region.get_startu(), region.get_endu()));
            for j in 0..bamlabels.len() {
                line.push_str(&format!("\t{}", matvec[[i, j]]));
            }
            writeln!(cfile, "{}", line).unwrap();
        }
    }
    if verbose {
        println!("Counts written.");
    }
    Ok(())
}

struct TempZip<I>
where I: Iterator {
    iterators: Vec<I>
}

impl<I, T> Iterator for TempZip<I>
where I: Iterator<Item=T> {
    type Item = Vec<T>;
    fn next(&mut self) -> Option<Self::Item> {
        let o: Option<Vec<T>> = self.iterators.iter_mut().map(|x| x.next()).collect();
        o
    }
}