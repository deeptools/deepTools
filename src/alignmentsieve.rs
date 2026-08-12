use crate::filehandler::{is_bed_or_gtf, read_bedfile};
use crate::filtering::Alignmentfilters;
use pyo3::prelude::*;
use pyo3::types::PyList;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::{self, Header, Read, Reader, Writer};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[pyfunction]
pub fn r_alignmentsieve(
    py: Python,
    bamifile: &str,               // input bamfile
    ofile: &str,                  // output file
    nproc: usize,                 // threads
    filter_metrics: &str,         // filter metrics file.
    filtered_out_readsfile: &str, // filtered_out_reads bam/bedfile.
    verbose: bool,                // verbose
    shift: Py<PyList>,            // python list of the shift to perform.
    bed: bool,                    // output format in BEDPE.
    filterrnastrand: &str,        // "forward", "reverse" or "None".
    minmappingquality: u8,        // minimum mapping quality.
    samflaginclude: u16,          // sam flag include
    samflagexclude: u16,          // sam flag exclude
    blacklist: &str,              // blacklist file name.
    minfraglen: u32,              // minimum fragment length.
    maxfraglen: u32,              // maximum fragment length.
    label: &str,                  // user defined label
    smartlabels: bool,            // derive label from filename
) -> PyResult<()> {
    // Open input BAM once and stream through in order
    let mut bam = Reader::from_path(bamifile).unwrap();
    // Reserve one extra thread for reader if there are enough threads
    let mut writerthreads = nproc;
    if nproc > 4 {
        let _ = bam.set_threads(2);
        writerthreads = nproc - 2;
    }
    let header = Header::from_template(bam.header());

    if verbose {
        println!("Opening BAM file: {}", bamifile);
    }

    let write_filters = filtered_out_readsfile != "None";
    let readshift: Vec<i32> = shift.extract(py).expect("Failed to extract shift");
    // shift is of length 0, 2, or 4.

    // Build chrom_names and chrom_sizes from BAM header — indexed by tid for O(1) lookup, zero allocations per record.
    let h = bam.header();
    let ntargets = h.target_count();
    let chrom_names: Vec<String> = (0..ntargets)
        .map(|tid| String::from_utf8(h.tid2name(tid).to_vec()).unwrap())
        .collect();
    let chrom_sizes: Vec<u64> = (0..ntargets)
        .map(|tid| h.target_len(tid).unwrap() as u64)
        .collect();

    // Load blacklist regions if provided
    let blacklist_regions = if blacklist != "None" {
        let isbed = is_bed_or_gtf(blacklist);
        let chrom_bounds: HashMap<String, u32> = chrom_names
            .iter()
            .cloned()
            .zip(chrom_sizes.iter().map(|&x| x as u32))
            .collect();
        match isbed.as_str() {
            "gtf" => panic!("Error: Please provide a bed file for the blacklist."),
            "bed" => {
                let (bls, _) = read_bedfile(&blacklist.to_string(), false, &chrom_bounds);
                Some(bls)
            }
            _ => panic!("Error: Cannot determine filetype of blacklist file."),
        }
    } else {
        None
    };

    let filters = Alignmentfilters::new(
        blacklist_regions,
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
        None,
    );

    // Open output writers
    let mut obam = if !bed {
        Some(Writer::from_path(ofile, &header, bam::Format::Bam).unwrap())
    } else {
        None
    };
    if let Some(ref mut w) = obam {
        let _ = w.set_threads(writerthreads);
    }
    let mut obed = if bed {
        Some(BufWriter::new(File::create(ofile).unwrap()))
    } else {
        None
    };

    let mut ofilterbam = if write_filters && !bed {
        Some(Writer::from_path(filtered_out_readsfile, &header, bam::Format::Bam).unwrap())
    } else {
        None
    };
    if let Some(ref mut w) = ofilterbam {
        let _ = w.set_threads(writerthreads);
    }
    let mut ofilterbed = if write_filters && bed {
        Some(BufWriter::new(
            File::create(filtered_out_readsfile).unwrap(),
        ))
    } else {
        None
    };

    let mut total_reads: u64 = 0;
    let mut filtered_reads: u64 = 0;

    // Lift the readshift branch outside the loop — eliminates a per-record branch
    // prediction target (readshift never changes during iteration).
    if !readshift.is_empty() {
        for result in bam.records() {
            let record = result.unwrap();
            total_reads += 1;

            let filtered = filter_record(&record, &chrom_names, &filters, record.tid() as usize);

            if filtered {
                filtered_reads += 1;
                if let Some(ref mut w) = ofilterbam {
                    w.write(&record).unwrap();
                } else if let Some(ref mut bw) = ofilterbed {
                    let tid = record.tid() as usize;
                    write_bed_line(&record, &chrom_names[tid], chrom_sizes[tid], bw);
                }
                continue;
            }

            let tid = record.tid() as usize;
            if let Some(shifted) = apply_shift(&record, &readshift, chrom_sizes[tid]) {
                if let Some(ref mut w) = obam {
                    w.write(&shifted).unwrap();
                } else if let Some(ref mut bw) = obed {
                    let tid = shifted.tid() as usize;
                    write_bed_line(&shifted, &chrom_names[tid], chrom_sizes[tid], bw);
                }
            }
        }
    } else {
        for result in bam.records() {
            let record = result.unwrap();
            total_reads += 1;

            let filtered = filter_record(&record, &chrom_names, &filters, record.tid() as usize);

            if filtered {
                filtered_reads += 1;
                if let Some(ref mut w) = ofilterbam {
                    w.write(&record).unwrap();
                } else if let Some(ref mut bw) = ofilterbed {
                    let tid = record.tid() as usize;
                    write_bed_line(&record, &chrom_names[tid], chrom_sizes[tid], bw);
                }
                continue;
            }

            if let Some(ref mut w) = obam {
                w.write(&record).unwrap();
            } else if let Some(ref mut bw) = obed {
                let tid = record.tid() as usize;
                write_bed_line(&record, &chrom_names[tid], chrom_sizes[tid], bw);
            }
        }
    }

    // Flush writers
    if let Some(mut w) = obed {
        w.flush().unwrap();
    }
    if let Some(mut w) = ofilterbed {
        w.flush().unwrap();
    }

    // Write filter metrics
    if filter_metrics != "None" {
        let sample_name = if smartlabels {
            smart_label(bamifile)
        } else if label != "None" {
            label.to_string()
        } else {
            bamifile.to_string()
        };

        let mut of = File::create(filter_metrics).unwrap();
        writeln!(of, "#bamFilterReads --filterMetrics").unwrap();
        writeln!(of, "#File\tReads Remaining\tTotal Initial Reads").unwrap();
        writeln!(
            of,
            "{}\t{}\t{}",
            sample_name,
            total_reads - filtered_reads,
            total_reads
        )
        .unwrap();
    }

    if verbose {
        println!(
            "Total reads: {}, Filtered: {}, Remaining: {}",
            total_reads,
            filtered_reads,
            total_reads - filtered_reads
        );
    }

    Ok(())
}

/// Determines whether a record should be filtered out.
/// true -> record should be filtered out.
fn filter_record(
    record: &bam::Record,
    chrom_names: &[String],
    alfilters: &Alignmentfilters,
    tid: usize,
) -> bool {
    let flags = record.flags();

    // Unmapped reads
    if flags & 4 != 0 {
        return true;
    }

    // Mapping quality
    if record.mapq() < alfilters.minmappingquality {
        return true;
    }

    // SAM flag include
    if alfilters.samflaginclude != 0
        && (flags & alfilters.samflaginclude) != alfilters.samflaginclude
    {
        return true;
    }

    // SAM flag exclude
    if alfilters.samflagexclude != 0 && (flags & alfilters.samflagexclude) != 0 {
        return true;
    }

    // Fragment length — check min and max independently
    // minFragmentLength > 0 and tLen < min; maxFragmentLength > 0 and tLen > max)
    if record.is_paired() {
        let tlen = record.insert_size().abs();
        if alfilters.minfraglen != 0 && tlen < alfilters.minfraglen as i64 {
            return true;
        }
        if alfilters.maxfraglen != 0 && tlen > alfilters.maxfraglen as i64 {
            return true;
        }
    } else if alfilters.minfraglen != 0 || alfilters.maxfraglen != 0 {
        let mut tlen: u32 = 0;
        for cig in record.cigar().iter() {
            match cig {
                bam::record::Cigar::Match(len)
                | bam::record::Cigar::Del(len)
                | bam::record::Cigar::Equal(len)
                | bam::record::Cigar::Diff(len) => tlen += len,
                _ => (),
            }
        }
        if alfilters.minfraglen != 0 && tlen < alfilters.minfraglen {
            return true;
        }
        if alfilters.maxfraglen != 0 && tlen > alfilters.maxfraglen {
            return true;
        }
    }

    // RNA strand filter
    if alfilters.filterrnastrand.as_str() != "None" {
        match (alfilters.filterrnastrand.as_str(), record.is_paired()) {
            ("forward", true) => {
                if !((flags & 144 == 128) || (flags & 96 == 64)) {
                    return true;
                }
            }
            ("forward", false) => {
                if !(flags & 16 == 16) {
                    return true;
                }
            }
            ("reverse", true) => {
                if !((flags & 144 == 144) || (flags & 96 == 96)) {
                    return true;
                }
            }
            ("reverse", false) => {
                if !(flags & 16 == 0) {
                    return true;
                }
            }
            _ => {}
        }
    }

    // Blacklist filter
    if alfilters.blacklist.is_some() {
        if alfilters.rec_in_blacklist(record, &chrom_names[tid]) {
            return true;
        }
    }

    false
}

/// Writes a single fragment as a BEDPE line (chrom\tstart\tend) if tlen > 0.
fn write_bed_line(record: &bam::Record, chrom: &str, chrom_len: u64, bw: &mut BufWriter<File>) {
    let tlen: i64 = {
        let ins = record.insert_size();
        if ins.abs() > 0 {
            ins
        } else {
            let mut cilen: i64 = 0;
            for cig in record.cigar().iter() {
                match cig {
                    bam::record::Cigar::Match(len)
                    | bam::record::Cigar::Del(len)
                    | bam::record::Cigar::Equal(len)
                    | bam::record::Cigar::Diff(len) => cilen += *len as i64,
                    _ => (),
                }
            }
            cilen
        }
    };

    // Only emit when tlen > 0 (keeps exactly one record per fragment for PE reads)
    if tlen <= 0 {
        return;
    }

    let start = record.pos() as i64;
    let mut end = start + tlen;

    if end > chrom_len as i64 {
        end = chrom_len as i64;
    }

    if end - start < 1 {
        return;
    }

    let _ = writeln!(bw, "{}\t{}\t{}", chrom, start, end);
}

// Calculate query_alignment_end: the 0-based end position of the alignment in the query,
// + excluding trailing soft-clips.
fn query_alignment_end(record: &bam::Record) -> i64 {
    let cigar: &[bam::record::Cigar] = &record.cigar().0;
    let mut read_pos: i64 = 0;
    for op in cigar {
        match op {
            bam::record::Cigar::Match(l)
            | bam::record::Cigar::Ins(l)
            | bam::record::Cigar::Equal(l)
            | bam::record::Cigar::Diff(l)
            | bam::record::Cigar::SoftClip(l) => read_pos += *l as i64,
            _ => (),
        }
    }
    // Remove trailing soft-clips
    for op in cigar.iter().rev() {
        if let bam::record::Cigar::SoftClip(l) = op {
            read_pos -= *l as i64;
        } else {
            break;
        }
    }
    read_pos
}

fn apply_shift(record: &bam::Record, shift: &[i32], chrom_len: u64) -> Option<bam::Record> {
    if !record.is_proper_pair() {
        return None;
    }
    if shift.len() < 4 {
        return None;
    }

    let tlen = record.insert_size();
    let mut start: i64 = record.pos();
    //  end = start + b.query_alignment_end
    let mut end: i64 = start + query_alignment_end(record);

    // is_first_in_template: true = read1, false = read2
    let is_read1 = record.is_first_in_template();
    let is_read2 = !is_read1;
    let is_rev = record.is_reverse();

    let delta_tlen: i64;
    if is_rev && is_read1 {
        end -= shift[2] as i64;
        delta_tlen = shift[3] as i64 - shift[2] as i64;
    } else if is_rev && is_read2 {
        end += shift[1] as i64;
        delta_tlen = shift[1] as i64 - shift[0] as i64;
    } else if !is_rev && is_read1 {
        start += shift[0] as i64;
        delta_tlen = shift[1] as i64 - shift[0] as i64;
    } else {
        start -= shift[3] as i64;
        delta_tlen = shift[3] as i64 - shift[2] as i64;
    }

    // Sanity check: if end - start < 1
    if end - start < 1 {
        if is_rev {
            start = end - 1;
        } else {
            end = start + 1;
        }
    }
    if start < 0 {
        start = 0;
    }
    if end > chrom_len as i64 {
        end = chrom_len as i64;
    }
    if end - start < 1 {
        return None;
    }

    let mut newrec = record.clone();

    let new_tlen = if tlen < 0 {
        tlen - delta_tlen
    } else {
        tlen + delta_tlen
    };

    // Build new CIGAR: single Match for fragment span ((0, end-start)).
    let span = (end - start) as u32;
    let new_cigar = CigarString(vec![bam::record::Cigar::Match(span)]);

    newrec.set(record.qname(), Some(&new_cigar), &[], &[]);

    // Remove all auxiliary tags
    let tags: Vec<Vec<u8>> = newrec
        .aux_iter()
        .filter_map(|r| r.ok().map(|(t, _)| t.to_vec()))
        .collect();
    for t in &tags {
        let _ = newrec.remove_aux(t);
    }

    newrec.set_pos(start);
    newrec.set_insert_size(new_tlen);

    // Mate position adjustment (Python: next_reference_start)
    let new_mpos = {
        let mut mpos = record.mpos();
        if is_read2 && is_rev {
            mpos += shift[0] as i64;
        } else if is_read1 && is_rev {
            mpos -= shift[3] as i64;
        }
        mpos
    };
    newrec.set_mpos(new_mpos);

    Some(newrec)
}

fn smart_label(label: &str) -> String {
    let basename = Path::new(label)
        .file_name()
        .and_then(|f| f.to_str())
        .unwrap_or(label);
    let without_ext = basename
        .split_once('.')
        .map(|(name, _)| name)
        .unwrap_or(basename);
    if without_ext.is_empty() {
        basename.to_string()
    } else {
        without_ext.to_string()
    }
}
