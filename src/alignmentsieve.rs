use crate::filehandler::{is_bed_or_gtf, read_bedfile};
use crate::filtering::Alignmentfilters;
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::PyList;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::{self, Header, Read, Reader, Writer};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[cfg(feature = "python")]
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
    let mut bam = Reader::from_path(bamifile)
        .unwrap_or_else(|e| panic!("Failed to open input BAM file '{}': {}", bamifile, e));
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
        .map(|tid| {
            String::from_utf8(h.tid2name(tid).to_vec()).unwrap_or_else(|e| {
                panic!(
                    "BAM header target name for tid {} in '{}' is not valid UTF-8: {}",
                    tid, bamifile, e
                )
            })
        })
        .collect();
    let chrom_sizes: Vec<u64> = (0..ntargets)
        .map(|tid| {
            h.target_len(tid).unwrap_or_else(|| {
                panic!(
                    "BAM header for '{}' is missing a target length for tid {}",
                    bamifile, tid
                )
            }) as u64
        })
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
        Some(
            Writer::from_path(ofile, &header, bam::Format::Bam)
                .unwrap_or_else(|e| panic!("Failed to create output BAM file '{}': {}", ofile, e)),
        )
    } else {
        None
    };
    if let Some(ref mut w) = obam {
        let _ = w.set_threads(writerthreads);
    }
    let mut obed = if bed {
        Some(BufWriter::new(File::create(ofile).unwrap_or_else(|e| {
            panic!("Failed to create output BED file '{}': {}", ofile, e)
        })))
    } else {
        None
    };

    let mut ofilterbam = if write_filters && !bed {
        Some(
            Writer::from_path(filtered_out_readsfile, &header, bam::Format::Bam).unwrap_or_else(
                |e| {
                    panic!(
                        "Failed to create filtered-out-reads BAM file '{}': {}",
                        filtered_out_readsfile, e
                    )
                },
            ),
        )
    } else {
        None
    };
    if let Some(ref mut w) = ofilterbam {
        let _ = w.set_threads(writerthreads);
    }
    let mut ofilterbed = if write_filters && bed {
        Some(BufWriter::new(
            File::create(filtered_out_readsfile).unwrap_or_else(|e| {
                panic!(
                    "Failed to create filtered-out-reads BED file '{}': {}",
                    filtered_out_readsfile, e
                )
            }),
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
            let record = result.unwrap_or_else(|e| {
                panic!(
                    "Failed to read a record from BAM file '{}': {}",
                    bamifile, e
                )
            });
            total_reads += 1;

            let filtered = filter_record(&record, &chrom_names, &filters, record.tid() as usize);

            if filtered {
                filtered_reads += 1;
                if let Some(ref mut w) = ofilterbam {
                    w.write(&record).unwrap_or_else(|e| {
                        panic!(
                            "Failed to write record '{}' to filtered-out-reads BAM file '{}': {}",
                            String::from_utf8_lossy(record.qname()),
                            filtered_out_readsfile,
                            e
                        )
                    });
                } else if let Some(ref mut bw) = ofilterbed {
                    let tid = record.tid() as usize;
                    write_bed_line(&record, &chrom_names[tid], chrom_sizes[tid], bw);
                }
                continue;
            }

            let tid = record.tid() as usize;
            if let Ok(shifted) = apply_shift(&record, &readshift, chrom_sizes[tid]) {
                if let Some(ref mut w) = obam {
                    w.write(&shifted).unwrap_or_else(|e| {
                        panic!(
                            "Failed to write shifted record '{}' to output BAM file '{}': {}",
                            String::from_utf8_lossy(shifted.qname()),
                            ofile,
                            e
                        )
                    });
                } else if let Some(ref mut bw) = obed {
                    let tid = shifted.tid() as usize;
                    write_bed_line(&shifted, &chrom_names[tid], chrom_sizes[tid], bw);
                }
            }
        }
    } else {
        for result in bam.records() {
            let record = result.unwrap_or_else(|e| {
                panic!(
                    "Failed to read a record from BAM file '{}': {}",
                    bamifile, e
                )
            });
            total_reads += 1;

            let filtered = filter_record(&record, &chrom_names, &filters, record.tid() as usize);

            if filtered {
                filtered_reads += 1;
                if let Some(ref mut w) = ofilterbam {
                    w.write(&record).unwrap_or_else(|e| {
                        panic!(
                            "Failed to write record '{}' to filtered-out-reads BAM file '{}': {}",
                            String::from_utf8_lossy(record.qname()),
                            filtered_out_readsfile,
                            e
                        )
                    });
                } else if let Some(ref mut bw) = ofilterbed {
                    let tid = record.tid() as usize;
                    write_bed_line(&record, &chrom_names[tid], chrom_sizes[tid], bw);
                }
                continue;
            }

            if let Some(ref mut w) = obam {
                w.write(&record).unwrap_or_else(|e| {
                    panic!(
                        "Failed to write record '{}' to output BAM file '{}': {}",
                        String::from_utf8_lossy(record.qname()),
                        ofile,
                        e
                    )
                });
            } else if let Some(ref mut bw) = obed {
                let tid = record.tid() as usize;
                write_bed_line(&record, &chrom_names[tid], chrom_sizes[tid], bw);
            }
        }
    }

    // Flush writers
    if let Some(mut w) = obed {
        w.flush()
            .unwrap_or_else(|e| panic!("Failed to flush output BED file '{}': {}", ofile, e));
    }
    if let Some(mut w) = ofilterbed {
        w.flush().unwrap_or_else(|e| {
            panic!(
                "Failed to flush filtered-out-reads BED file '{}': {}",
                filtered_out_readsfile, e
            )
        });
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

        let mut of = File::create(filter_metrics).unwrap_or_else(|e| {
            panic!(
                "Failed to create filter metrics file '{}': {}",
                filter_metrics, e
            )
        });
        writeln!(of, "#bamFilterReads --filterMetrics").unwrap_or_else(|e| {
            panic!(
                "Failed to write to filter metrics file '{}': {}",
                filter_metrics, e
            )
        });
        writeln!(of, "#File\tReads Remaining\tTotal Initial Reads").unwrap_or_else(|e| {
            panic!(
                "Failed to write to filter metrics file '{}': {}",
                filter_metrics, e
            )
        });
        writeln!(
            of,
            "{}\t{}\t{}",
            sample_name,
            total_reads - filtered_reads,
            total_reads
        )
        .unwrap_or_else(|e| {
            panic!(
                "Failed to write to filter metrics file '{}': {}",
                filter_metrics, e
            )
        });
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

fn cigar_op_len(op: &bam::record::Cigar) -> u32 {
    use bam::record::Cigar::*;
    match *op {
        Match(n) | Ins(n) | Del(n) | RefSkip(n) | SoftClip(n) | HardClip(n) | Pad(n) | Equal(n) | Diff(n) => n,
    }
}

fn cigar_with_len(op: bam::record::Cigar, n: u32) -> bam::record::Cigar {
    use bam::record::Cigar::*;
    match op {
        Match(_) => Match(n),
        Ins(_) => Ins(n),
        Del(_) => Del(n),
        RefSkip(_) => RefSkip(n),
        SoftClip(_) => SoftClip(n),
        HardClip(_) => HardClip(n),
        Pad(_) => Pad(n),
        Equal(_) => Equal(n),
        Diff(_) => Diff(n),
    }
}

fn cigar_consumes_ref(op: &bam::record::Cigar) -> bool {
    use bam::record::Cigar::*;
    matches!(op, Match(_) | Equal(_) | Diff(_) | Del(_) | RefSkip(_))
}

/// Remove `k` reference bases from the left of a CIGAR. Query bases removed from
/// M/=/X (and any I/S before them) become one leading soft clip, so the query-consuming
/// length, and therefore SEQ/QUAL, is unchanged.
fn trim_cigar_left(ops: &[bam::record::Cigar], mut k: u32) -> Vec<bam::record::Cigar> {
    use bam::record::Cigar::*;
    let mut q = 0u32;
    let mut rest: Vec<bam::record::Cigar> = Vec::new();
    for (i, op) in ops.iter().enumerate() {
        if k == 0 {
            rest.extend_from_slice(&ops[i..]);
            break;
        }
        let n = cigar_op_len(op);
        if !cigar_consumes_ref(op) {
            if matches!(op, Ins(_) | SoftClip(_)) {
                q += n;
            }
            continue;
        }
        let t = n.min(k);
        if matches!(op, Match(_) | Equal(_) | Diff(_)) {
            q += t;
        }
        k -= t;
        if t < n {
            rest.push(cigar_with_len(*op, n - t));
            rest.extend_from_slice(&ops[i + 1..]);
            break;
        }
    }
    let mut out = Vec::with_capacity(rest.len() + 1);
    if q > 0 {
        out.push(SoftClip(q));
    }
    out.extend(rest);
    out
}

/// Why `apply_shift` refused a record. Lets callers count each reason separately.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShiftRefusal {
    /// Record is not part of a proper pair.
    NotProperPair,
    /// Fewer than 4 shift values were given.
    ShortShift,
    /// Shift leaves less than 1 bp of span after clamping to the chromosome.
    ZeroSpan,
    /// Shift would move an alignment end past its original span (start left or end right).
    LeftExtension,
}

pub fn apply_shift(record: &bam::Record, shift: &[i32], chrom_len: u64) -> Result<bam::Record, ShiftRefusal> {
    if !record.is_proper_pair() {
        return Err(ShiftRefusal::NotProperPair);
    }
    if shift.len() < 4 {
        return Err(ShiftRefusal::ShortShift);
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
        return Err(ShiftRefusal::ZeroSpan);
    }

    let mut newrec = record.clone();

    let new_tlen = if tlen < 0 {
        tlen - delta_tlen
    } else {
        tlen + delta_tlen
    };

    // Keep SEQ, QUAL and aux tags. Only the alignment moves: trim reference bases
    // from the left (shift start right) and/or the right (shift end left).
    // Extending a read would need reference bases we do not have, so it is refused.
    let old_start = record.pos();
    let old_end = old_start + query_alignment_end(record);
    let d_start = start - old_start;
    let d_end = end - old_end;
    if d_start < 0 || d_end > 0 {
        return Err(ShiftRefusal::LeftExtension);
    }
    let ops: Vec<bam::record::Cigar> = record.cigar().iter().cloned().collect();
    let left = trim_cigar_left(&ops, d_start as u32);
    let mut rev: Vec<bam::record::Cigar> = left.into_iter().rev().collect();
    rev = trim_cigar_left(&rev, (-d_end) as u32);
    let new_ops: Vec<bam::record::Cigar> = rev.into_iter().rev().collect();
    let seq: Vec<u8> = record.seq().as_bytes();
    let qual: Vec<u8> = record.qual().to_vec();
    newrec.set(record.qname(), Some(&CigarString(new_ops)), &seq, &qual);

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

    Ok(newrec)
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

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Aux, Cigar};

    fn proper_pair(cigar: &[Cigar], pos: i64) -> bam::Record {
        let qlen: u32 = cigar
            .iter()
            .map(|c| match c {
                Cigar::Match(n) | Cigar::SoftClip(n) | Cigar::Ins(n) | Cigar::Equal(n) | Cigar::Diff(n) => *n,
                _ => 0,
            })
            .sum();
        let seq = vec![b'A'; qlen as usize];
        let qual = vec![30u8; qlen as usize];
        let mut r = bam::Record::new();
        r.set(b"read1", Some(&CigarString(cigar.to_vec())), &seq, &qual);
        r.set_pos(pos);
        r.set_tid(0);
        r.set_mtid(0);
        r.set_mpos(pos + 200);
        r.set_insert_size(300);
        r.set_flags(0x1 | 0x2 | 0x40);
        r.push_aux(b"CB", Aux::String("AAAC-1")).unwrap();
        r
    }

    #[test]
    fn shift_keeps_cell_barcode_tag() {
        let r = proper_pair(&[Cigar::Match(50)], 1000);
        let out = apply_shift(&r, &[4, -5, 5, -4], 1_000_000).expect("shifted");
        assert!(out.aux(b"CB").is_ok(), "CB tag dropped by apply_shift");
    }

    #[test]
    fn shift_keeps_soft_clip_in_cigar() {
        let r = proper_pair(&[Cigar::SoftClip(5), Cigar::Match(50)], 1000);
        let out = apply_shift(&r, &[4, -5, 5, -4], 1_000_000).expect("shifted");
        assert_eq!(out.pos(), 1004);
        assert_eq!(out.cigar().end_pos() - out.pos(), 46, "ref span must equal end - start");
        let q: u32 = out.cigar().iter().filter(|c| !matches!(c, Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::HardClip(_) | Cigar::Pad(_))).map(|c| match c { Cigar::Match(n) | Cigar::Ins(n) | Cigar::SoftClip(n) | Cigar::Equal(n) | Cigar::Diff(n) => *n, _ => 0 }).sum();
        assert_eq!(q as usize, out.seq_len(), "CIGAR query length must equal SEQ length");
    }
}

#[cfg(test)]
mod refusal_tests {
    use super::*;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};

    fn fwd_read1(pos: i64) -> bam::Record {
        let seq = vec![b'A'; 50];
        let qual = vec![30u8; 50];
        let mut r = bam::Record::new();
        r.set(b"read1", Some(&CigarString(vec![Cigar::Match(50)])), &seq, &qual);
        r.set_pos(pos);
        r.set_tid(0);
        r.set_mtid(0);
        r.set_mpos(pos + 200);
        r.set_insert_size(300);
        r.set_flags(0x1 | 0x2 | 0x40);
        r.push_aux(b"CB", Aux::String("AAAC-1")).unwrap();
        r
    }

    #[test]
    fn forward_read1_negative_start_shift_is_left_extension() {
        let r = fwd_read1(1000);
        assert_eq!(
            apply_shift(&r, &[-3, 0, 0, 0], 1_000_000).unwrap_err(),
            ShiftRefusal::LeftExtension
        );
    }
}
