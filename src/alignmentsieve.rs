use crate::filehandler::{is_bed_or_gtf, read_bedfile};
use crate::filtering::Alignmentfilters;
use pyo3::prelude::*;
use pyo3::types::PyList;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::{self, Header, HeaderView, Read, Reader, Writer};
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
    let header = Header::from_template(bam.header());
    let header_view = bam.header().clone();

    if verbose {
        println!("Opening BAM file: {}", bamifile);
    }

    let write_filters = filtered_out_readsfile != "None";
    let readshift: Vec<i32> = shift.extract(py).expect("Failed to extract shift");
    // shift is of length 0, 2, or 4.

    // Build chrom_sizes from BAM header
    let chrom_sizes: HashMap<String, u64> = (0..header_view.target_count())
        .map(|tid| {
            (
                String::from_utf8(header_view.tid2name(tid).to_vec()).unwrap(),
                header_view.target_len(tid).unwrap() as u64,
            )
        })
        .collect();

    // Load blacklist regions if provided
    let blacklist_regions = if blacklist != "None" {
        let isbed = is_bed_or_gtf(blacklist);
        let chrom_keys: Vec<&String> = chrom_sizes.keys().collect();
        match isbed.as_str() {
            "gtf" => panic!("Error: Please provide a bed file for the blacklist."),
            "bed" => {
                let (bls, _) = read_bedfile(&blacklist.to_string(), false, chrom_keys);
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
        let _ = w.set_threads(nproc);
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
        let _ = w.set_threads(nproc);
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

    for result in bam.records() {
        let record = result.unwrap();
        total_reads += 1;

        let filtered = filter_record(&record, &header_view, &filters);

        if filtered {
            filtered_reads += 1;
            if let Some(ref mut w) = ofilterbam {
                w.write(&record).unwrap();
            } else if let Some(ref mut bw) = ofilterbed {
                write_bed_line(&record, &header_view, &chrom_sizes, bw);
            }
            continue;
        }

        if !readshift.is_empty() {
            let chrom_name =
                String::from_utf8(header_view.tid2name(record.tid() as u32).to_vec()).unwrap();
            if let Some(shifted) =
                apply_shift(&record, &readshift, chrom_sizes.get(&chrom_name).copied())
            {
                if let Some(ref mut w) = obam {
                    w.write(&shifted).unwrap();
                } else if let Some(ref mut bw) = obed {
                    write_bed_line(&shifted, &header_view, &chrom_sizes, bw);
                }
            }
            continue;
        }

        // No shift — write original record
        if let Some(ref mut w) = obam {
            w.write(&record).unwrap();
        } else if let Some(ref mut bw) = obed {
            write_bed_line(&record, &header_view, &chrom_sizes, bw);
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
/// Returns true if the record fails any filter criterion.
fn filter_record(
    record: &bam::Record,
    bam_header: &HeaderView,
    alfilters: &Alignmentfilters,
) -> bool {
    // Unmapped reads
    if record.is_unmapped() {
        return true;
    }

    // Mapping quality
    if record.mapq() < alfilters.minmappingquality {
        return true;
    }

    // SAM flag include
    if alfilters.samflaginclude != 0
        && (record.flags() & alfilters.samflaginclude) != alfilters.samflaginclude
    {
        return true;
    }

    // SAM flag exclude
    if alfilters.samflagexclude != 0 && (record.flags() & alfilters.samflagexclude) != 0 {
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
    } else {
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
                if !((record.flags() & 144 == 128) || (record.flags() & 96 == 64)) {
                    return true;
                }
            }
            ("forward", false) => {
                if !(record.flags() & 16 == 16) {
                    return true;
                }
            }
            ("reverse", true) => {
                if !((record.flags() & 144 == 144) || (record.flags() & 96 == 96)) {
                    return true;
                }
            }
            ("reverse", false) => {
                if !(record.flags() & 16 == 0) {
                    return true;
                }
            }
            _ => {}
        }
    }

    // Blacklist filter
    if alfilters.blacklist.is_some() {
        let chrom_name =
            String::from_utf8(bam_header.tid2name(record.tid() as u32).to_vec()).unwrap();
        if alfilters.rec_in_blacklist(record, &chrom_name) {
            return true;
        }
    }

    false
}

/// Writes a single fragment as a BEDPE line (chrom\tstart\tend) if tlen > 0.
fn write_bed_line(
    record: &bam::Record,
    bam_header: &HeaderView,
    chrom_sizes: &HashMap<String, u64>,
    bw: &mut BufWriter<File>,
) {
    let tlen: i64 = if record.insert_size().abs() > 0 {
        record.insert_size()
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
    };

    // Only emit when tlen > 0 (keeps exactly one record per fragment for PE reads)
    if tlen <= 0 {
        return;
    }

    let chrom = String::from_utf8(bam_header.tid2name(record.tid() as u32).to_vec()).unwrap();
    let start = record.pos() as i64;
    let mut end = start + tlen;

    if let Some(&chrom_len) = chrom_sizes.get(&chrom) {
        if end > chrom_len as i64 {
            end = chrom_len as i64;
        }
    }

    if end - start < 1 {
        return;
    }

    let _ = writeln!(bw, "{}\t{}\t{}", chrom, start, end);
}

// Calculate query_alignment_end: the 0-based end position of the alignment in the query,
// + excluding trailing soft-clips.
fn query_alignment_end(record: &bam::Record) -> i64 {
    let cigar_ops: Vec<bam::record::Cigar> = record.cigar().iter().cloned().collect();
    let mut read_pos: i64 = 0;
    for op in &cigar_ops {
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
    for op in cigar_ops.iter().rev() {
        if let bam::record::Cigar::SoftClip(l) = op {
            read_pos -= *l as i64;
        } else {
            break;
        }
    }
    read_pos
}

fn apply_shift(record: &bam::Record, shift: &[i32], chrom_len: Option<u64>) -> Option<bam::Record> {
    if !record.is_proper_pair() {
        return None;
    }
    if shift.len() < 4 {
        return None;
    }

    let tLen = record.insert_size();
    let mut start: i64 = record.pos();
    //  end = start + b.query_alignment_end
    let mut end: i64 = start + query_alignment_end(record);

    // is_first_in_template: true = read1, false = read2
    let is_read1 = record.is_first_in_template();
    let is_read2 = !is_read1;
    let is_rev = record.is_reverse();

    let mut deltaTLen: i64 = 0;
    if is_rev && is_read1 {
        end -= shift[2] as i64;
        deltaTLen = shift[3] as i64 - shift[2] as i64;
    } else if is_rev && is_read2 {
        end += shift[1] as i64;
        deltaTLen = shift[1] as i64 - shift[0] as i64;
    } else if !is_rev && is_read1 {
        start += shift[0] as i64;
        deltaTLen = shift[1] as i64 - shift[0] as i64;
    } else {
        start -= shift[3] as i64;
        deltaTLen = shift[3] as i64 - shift[2] as i64;
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
    if let Some(len) = chrom_len {
        if end > len as i64 {
            end = len as i64;
        }
    }
    if end - start < 1 {
        return None;
    }

    let mut newrec = record.clone();

    // Template length: tLen + deltaTLen (or tLen - deltaTLen if tLen < 0)
    let new_tlen = if tLen < 0 {
        tLen - deltaTLen
    } else {
        tLen + deltaTLen
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
