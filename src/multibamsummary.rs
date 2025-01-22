use pyo3::prelude::*;
use pyo3::types::PyList;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use itertools::{multizip, multiunzip, izip};
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use ndarray::{Array1, Array2, stack};
use std::io;
use ndarray_npy::NpzWriter;
use std::borrow::Cow;
use crate::covcalc::{bam_pileup, parse_regions, Alignmentfilters};
use crate::filehandler::{bam_ispaired};
use crate::calc::{median, calc_ratio};
use crate::bamcompare::ParsedBamFile;
use crate::normalization::scale_factor_bamcompare;

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
    binsize: u32,
    distance_between_bins: u32,
    nproc: usize,
    _bed_file: &str,
    regions: Vec<(String, u32, u32)>,
    _blacklist: &str,
    _verbose: bool,
    _extend_reads: u32,
    _center_reads: bool,
    sam_flag_incl: u16, // sam flag include
    sam_flag_excl: u16, // sam flag exclude
    min_fragment_length: u32, // minimum fragment length.
    max_fragment_length: u32, // maximum fragment length.
    min_mapping_quality: u8, // minimum mapping quality.
    _keep_exons: bool,
    _txnid: &str, // transcript id to use when parsing GTF file
    _exonid: &str, // exon id to use when parsing GTF file
    _txniddesignator: &str, // designator to use when parsing GTF file
) -> PyResult<()> {
    let mut bamfiles: Vec<String> = Vec::new();
    let mut bamlabels: Vec<String> = Vec::new();
    let mut ignorechr: Vec<String> = Vec::new();
    Python::with_gil(|py| {
        bamfiles = bam_files.extract(py).expect("Failed to retrieve bam files.");
        bamlabels = labels.extract(py).expect("Failed to retrieve labels.");
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
    for (_ispe, _bf) in ispe.iter().zip(bamfiles.iter()) {
        println!("Sample: {} is-paired: {}", _bf, _ispe);
    }

    let filters: Alignmentfilters = Alignmentfilters {
        minmappingquality: min_mapping_quality,
        samflaginclude: sam_flag_incl,
        samflagexclude: sam_flag_excl,
        minfraglen: min_fragment_length,
        maxfraglen: max_fragment_length
    };
    
    let (regions, chromsizes) = parse_regions(&regions, bamfiles.get(0).unwrap());
    let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();

    
    // Zip together bamfiles and ispe into a vec of tuples.
    let bampfiles: Vec<_> = bamfiles.into_iter().zip(ispe.into_iter()).collect();

    let covcalcs: Vec<_> = pool.install(|| {
        bampfiles.par_iter()
            .map(|(bamfile, ispe)| {
                let (bg, mapped, unmapped, readlen, fraglen) = regions.par_iter()
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

    
    // Scale factors (when needed)
    // relevant python code. 
    // loggeomeans = np.sum(np.log(m), axis=1) / m.shape[1]
    // # Mask after computing the geometric mean
    // m = np.ma.masked_where(m <= 0, m)
    // loggeomeans = np.ma.masked_where(np.isinf(loggeomeans), loggeomeans)
    // # DESeq2 ratio-based size factor
    // sf = np.exp(np.ma.median((np.log(m).T - loggeomeans).T, axis=0))
    // return 1. / sf
    // create output file to write to

    // create output file to write lines into
    let mut cfile = None;
    if out_raw_counts != "None" {
        cfile = Some(File::create(out_raw_counts).unwrap());
    }
    
    let mut matvec: Vec<Array1<f32>> = Vec::new();
    
    let its: Vec<_> = covcalcs.iter().map(|x| x.iter()).collect();
    let zips = TempZip { iterators: its };
    for c in zips {
        // set up Readers
        let readers: Vec<_> = c.iter().map(|x| BufReader::new(File::open(x).unwrap()).lines()).collect();
        
        for mut _l in (TempZip { iterators: readers }) {
            // unwrap all lines in _l
            let lines: Vec<_> = _l
                .iter_mut()
                .map(|x| x.as_mut().unwrap())
                .map(|x| x.split('\t').collect())
                .map(|x: Vec<&str> | ( x[0].to_string(), x[1].parse::<u32>().unwrap(), x[2].parse::<u32>().unwrap(), x[3].parse::<f32>().unwrap() ) )
                .collect();
            let arrline = Array1::from(lines.iter().map(|x| x.3).collect::<Vec<_>>());
            matvec.push(arrline);
            // Write to countsfile if cfile is not None
            if let Some(ref mut cfile) = cfile {
                let mut line: Vec<Countline> = Vec::new();
                line.push(Countline::Text(lines[0].0.clone()));
                line.push(Countline::Int(lines[0].1));
                line.push(Countline::Int(lines[0].2));
                for _line in lines {
                    line.push(Countline::Float(_line.3));
                }
                let fline = line.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\t");
                writeln!(cfile, "{}", fline).unwrap();
            }
        }
    }
    let array2: Array2<f32> = Array2::from_shape_vec(
        (matvec.len(), matvec[0].len()),
        matvec.into_iter().flat_map(|x| x.into_iter()).collect()
    ).unwrap();
    let mut npz = NpzWriter::new(File::create(ofile).unwrap());
    npz.add_array("matrix", &array2).unwrap();
    npz.add_array("labels", &bamlabels_arr).unwrap();
    npz.finish().unwrap();

    Ok(())
}

#[derive(Debug)]
enum Countline {
    Int(u32),
    Float(f32),
    Text(String),
}

impl Countline {
    fn to_string(&self) -> String {
        match self {
            Countline::Int(i) => i.to_string(),
            Countline::Float(f) => f.to_string(),
            Countline::Text(t) => t.clone(),
        }
    }
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