use pyo3::prelude::*;
use pyo3::types::PyList;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use rust_htslib::bam::{self, Header, IndexedReader, Read, Reader, Writer};
use tempfile::{Builder, TempPath};
use std::fs::File;
use std::io::Write;
use crate::covcalc::{parse_regions, Region};
use crate::filtering::Alignmentfilters;
use crate::filehandler::{is_bed_or_gtf, read_bedfile};

#[pyfunction]
pub fn r_alignmentsieve(
    bamifile: &str, // input bamfile
    ofile: &str, // output file
    nproc: usize, // threads
    filter_metrics: &str, // filter metrics file.
    filtered_out_readsfile: &str, // filtered_out_reads bam/bedfile.
    verbose: bool, // verbose
    shift: Py<PyList>, // python list of the shift to perform.
    _bed: bool, // output format in BEDPE.
    filterrnastrand: &str, // "forward", "reverse" or "None".
    minmappingquality: u8, // minimum mapping quality.
    samflaginclude: u16, // sam flag include
    samflagexclude: u16, // sam flag exclude
    blacklist: &str, // blacklist file name.
    minfraglen: u32, // minimum fragment length.
    maxfraglen: u32, // maximum fragment length.
    _extend_reads: u32,
    _center_reads: bool,
    
) -> PyResult<()> {
    // Input bam file
    let bam = Reader::from_path(bamifile).unwrap();
    let header = Header::from_template(bam.header());
    let _header_view = bam.header().clone();

    let mut write_filters: bool = false;
    if filtered_out_readsfile != "None" {
        write_filters = true;
    }
    let mut readshift: Vec<i32> = Vec::new();
    Python::with_gil(|py| {
        readshift = shift.extract(py).expect("Failed to extract shift.");
    });
    // shift is of length 0, 2, or 4.

    // Define regions 
    let (regions, chromsizes) = parse_regions("None", vec![bamifile]);
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

    let filters = Alignmentfilters::new(
        backlistregions,
        Some(minmappingquality),
        Some(samflaginclude),
        Some(samflagexclude),
        Some(minfraglen),
        Some(maxfraglen),
        None,
        None,
        Some(filterrnastrand.to_string()),
        None,
        None,
    );
    let pool = ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    let (sieve, filtersieve, totalreads, filteredreads) = pool.install(|| {
        regions.par_iter()
            .map(|i| sieve_bamregion(bamifile, i, &filters, &readshift, write_filters, nproc, verbose))
            .reduce(
                || (Vec::new(), Vec::new(), 0, 0),
                |(mut _sieve, mut _filtersieve, mut _total, mut _filter), (sieve, filtersieve, total, filter)| {
                    _sieve.extend(sieve);
                    _filtersieve.extend(filtersieve);
                    _total += total;
                    _filter += filter;
                    (_sieve, _filtersieve, _total, _filter)
                }
            )
    });

    // write output
    let mut obam = Writer::from_path(ofile, &header, bam::Format::Bam).unwrap();
    let _ = obam.set_threads(nproc);
    for sb in sieve.into_iter() {
        if let Some(sb) = sb {
            let mut bam = Reader::from_path(&sb).unwrap();
            for result in bam.records() {
                let record = result.unwrap();
                obam.write(&record).unwrap();
            }
        }
    }
    // write filtered reads if necessary
    if write_filters {
        let mut ofilterbam = Writer::from_path(filtered_out_readsfile, &header, bam::Format::Bam).unwrap();
        let _ = ofilterbam.set_threads(nproc);
        for sb in filtersieve.into_iter() {
            if let Some(sb) = sb {
                let mut bam = Reader::from_path(&sb).unwrap();
                for result in bam.records() {
                    let record = result.unwrap();
                    ofilterbam.write(&record).unwrap();
                }
            }
        }
    }

    let _ofilterbam = Writer::from_path(filtered_out_readsfile, &header, bam::Format::Bam).unwrap();

    if filter_metrics != "None" {
        let mut of = File::create(filter_metrics).unwrap();
        // write header
        writeln!(of, "#bamFilterReads --filterMetrics").unwrap();
        writeln!(of, "#File\tReads\tRemaining Total\tInitial Reads").unwrap();
        writeln!(of, "{}\t{}\t{}", bamifile, totalreads-filteredreads, totalreads).unwrap();
    }

    Ok(())
}


fn sieve_bamregion(ibam: &str, regstruct: &Region, alfilters: &Alignmentfilters, _shift: &Vec<i32>, write_filters: bool, nproc: usize, verbose: bool) -> (Vec<Option<TempPath>>, Vec<Option<TempPath>>, u64, u64) {
    let region = (regstruct.chrom.clone(), regstruct.get_startu(), regstruct.get_endu());
    let mut total_reads: u64 = 0;
    let mut filtered_reads: u64 = 0;
    let mut bam = IndexedReader::from_path(ibam).unwrap();
    let header = Header::from_template(bam.header());

    let mut written = false;
    let mut filterwritten = false;

    let sievebam = Builder::new()
        .prefix("deeptoolstmp_alsieve_")
        .suffix(".bam")
        .rand_bytes(12)
        .tempfile()
        .expect("Failed to create temporary file.");

    let sievebam_path = sievebam.into_temp_path();
    let mut sievebamout = Writer::from_path(&sievebam_path, &header, bam::Format::Bam).unwrap();

    let filterbam = Builder::new()
        .prefix("deeptoolstmp_alsieve_filtered_")
        .suffix(".bam")
        .rand_bytes(12)
        .tempfile()
        .expect("Failed to create temporary file.");
    let filterbam_path = filterbam.into_temp_path();
    let mut filterbamout = if write_filters {
        Some(Writer::from_path(&filterbam_path, &header, bam::Format::Bam).unwrap())
    } else {
        None
    };
    if nproc > 4 {
        let readthreads = 2;
        let writethreads = nproc - 2;
        let _ = bam.set_threads(readthreads);
        let _ = sievebamout.set_threads(writethreads);
        if verbose {
            println!("Reading = {}, Writing = {}", readthreads, writethreads);
        }
    }

    bam.fetch((region.0.as_str(), region.1, region.2)).unwrap();

    for result in bam.records() {
        let record = result.unwrap();
        total_reads += 1;

        // Filter reads
        // Filter unmapped reads.
        if record.is_unmapped() {
            filtered_reads += 1;
            if let Some(filterbamout) = &mut filterbamout {
                filterbamout.write(&record).unwrap();
                filterwritten = true;
            }
            continue;
        }
        // Mapping qualities.
        if record.mapq() < alfilters.minmappingquality {
            filtered_reads += 1;
            if let Some(filterbamout) = &mut filterbamout {
                filterbamout.write(&record).unwrap();
                filterwritten = true;
            }
            continue;
        }
        
        // SAM flags
        if alfilters.samflaginclude != 0 && (record.flags() & alfilters.samflaginclude) == 0 {
            filtered_reads += 1;
            if let Some(filterbamout) = &mut filterbamout {
                filterbamout.write(&record).unwrap();
                filterwritten = true;
            }
            continue;
        }
        if alfilters.samflagexclude != 0 && (record.flags() & alfilters.samflagexclude) != 0 {
            filtered_reads += 1;
            if let Some(filterbamout) = &mut filterbamout {
                filterbamout.write(&record).unwrap();
                filterwritten = true;
            }
            continue;
        }

        // fragment length
        if alfilters.minfraglen != 0 || alfilters.maxfraglen != 0 {
            if record.is_paired() {
                if record.insert_size().abs() < alfilters.minfraglen as i64 || record.insert_size().abs() > alfilters.maxfraglen as i64 {
                    filtered_reads += 1;
                    if let Some(filterbamout) = &mut filterbamout {
                        filterbamout.write(&record).unwrap();
                        filterwritten = true;
                    }
                    continue;
                }
            } else {
                // Parse cigartuples
                let mut tlen: u32 = 0;
                for cig in record.cigar().iter() {
                    match cig {
                        bam::record::Cigar::Match(len) => tlen += len,
                        bam::record::Cigar::Del(len) => tlen += len,
                        bam::record::Cigar::Equal(len) => tlen += len,
                        bam::record::Cigar::Diff(len) => tlen += len,
                        _ => (),
                    }
                }
                if tlen < alfilters.minfraglen || tlen > alfilters.maxfraglen {
                    filtered_reads += 1;
                    if let Some(filterbamout) = &mut filterbamout {
                        filterbamout.write(&record).unwrap();
                        filterwritten = true;
                    }
                    continue;
                }
            }
        }
        if alfilters.filterrnastrand.as_str() != "None" {
            match (alfilters.filterrnastrand.as_str(), record.is_paired()) {
                ("forward", true) => {
                    if !((record.flags() & 144 == 128) || (record.flags() & 96 == 64)) {
                        filtered_reads += 1;
                        if let Some(filterbamout) = &mut filterbamout {
                            filterbamout.write(&record).unwrap();
                            filterwritten = true;
                        }
                        continue;
                    }
                },
                ("forward", false) => {
                    if !(record.flags() & 16 == 16) {
                        filtered_reads += 1;
                        if let Some(filterbamout) = &mut filterbamout {
                            filterbamout.write(&record).unwrap();
                            filterwritten = true;
                        }
                        continue;
                    }
                },
                ("reverse", true) => {
                    if !((record.flags() & 144 == 144) || (record.flags() & 96 == 96)) {
                        filtered_reads += 1;
                        if let Some(filterbamout) = &mut filterbamout {
                            filterbamout.write(&record).unwrap();
                            filterwritten = true;
                        }
                        continue;
                    }
                },
                ("reverse", false) => {
                    if !(record.flags() & 16 == 0) {
                        filtered_reads += 1;
                        if let Some(filterbamout) = &mut filterbamout {
                            filterbamout.write(&record).unwrap();
                            filterwritten = true;
                        }
                        continue;
                    }
                },
                _ => {},
            }
        }
        sievebamout.write(&record).unwrap();
        written = true;
    }

    match (written, filterwritten) {
        (true, true) => {
            (vec![Some(sievebam_path)], vec![Some(filterbam_path)], total_reads, filtered_reads)
        },
        (true, false) => {
            (vec![Some(sievebam_path)], vec![None], total_reads, filtered_reads)
        },
        (false, true) => {
            (vec![None], vec![Some(filterbam_path)], total_reads, filtered_reads)
        },
        (false, false) => {
            (vec![None], vec![None], total_reads, filtered_reads)
        }
    }

}