use crate::calc::median;
use crate::covcalc::Region;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use rust_htslib::bam::{IndexedReader, Read, Record, ext::BamRecordExtensions};
use std::collections::HashMap;

#[derive(Clone)]
pub struct BlacklistIndex {
    // Maps chromosome name to sorted list of (start, end) tuples
    chroms: HashMap<String, Vec<(u32, u32)>>,
}

impl BlacklistIndex {
    pub fn from_regions(regions: &[Region]) -> Self {
        let mut chroms: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
        for r in regions {
            chroms
                .entry(r.chrom.clone())
                .or_default()
                .push((r.get_startu(), r.get_endu()));
        }
        for v in chroms.values_mut() {
            v.sort_unstable_by_key(|&(s, _)| s);
        }
        Self { chroms }
    }

    pub fn contains(&self, chrom: &str, pos: u32) -> bool {
        let list = match self.chroms.get(chrom) {
            Some(l) => l,
            None => return false,
        };
        match list.binary_search_by_key(&pos, |&(s, _)| s) {
            Ok(ix) => list[ix].1 >= pos,
            Err(ix) => {
                if ix > 0 && list[ix - 1].0 <= pos && list[ix - 1].1 >= pos {
                    return true;
                }
                if ix < list.len() && list[ix].0 <= pos && list[ix].1 >= pos {
                    return true;
                }
                false
            }
        }
    }

    // pub fn overlaps(&self, chrom: &str, start: u32, end: u32) -> bool {
    //     // Returns true if any position in [start, end) falls within a blacklisted region.
    //     let list = match self.chroms.get(chrom) {
    //         Some(l) => l,
    //         None => return false,
    //     };
    //     // Find the first blacklist entry with start >= start (or the one just before it)
    //     match list.binary_search_by_key(&start, |&(s, _)| s) {
    //         Ok(ix) => list[ix].1 > start,
    //         Err(ix) => {
    //             // Check the entry just before ix (its start < `start`), see if it extends past `start`
    //             if ix > 0 && list[ix - 1].1 > start {
    //                 return true;
    //             }
    //             // Check entries starting from ix onward, as long as their start < end
    //             if ix < list.len() && list[ix].0 < end {
    //                 return true;
    //             }
    //             false
    //         }
    //     }
    // }
}

#[derive(Clone)]
pub struct Alignmentfilters {
    pub blacklist: Option<Vec<Region>>,
    blacklist_index: Option<BlacklistIndex>,
    pub minmappingquality: u8,
    pub samflaginclude: u16,
    pub samflagexclude: u16,
    pub minfraglen: u32,
    pub maxfraglen: u32,
    pub mnase: bool,
    pub offset: (i32, i32),
    pub filterrnastrand: String,
    pub extendreads: bool,
    pub extendreadslen: u32,
    pub extendreads_auto: bool,
    pub centerreads: bool,
    pub filter: bool,
    pub manipulate: bool,
}
impl Alignmentfilters {
    pub fn new(
        blacklist: Option<Vec<Region>>,
        minmappingquality: Option<u8>,
        samflaginclude: Option<u16>,
        samflagexclude: Option<u16>,
        minfraglen: Option<u32>,
        maxfraglen: Option<u32>,
        mnase: Option<bool>,
        offset: Option<(i32, i32)>,
        filterrnastrand: Option<String>,
        extendreads: Option<bool>,
        extendreadslen: Option<u32>,
        centerreads: Option<bool>,
    ) -> Self {
        // Go through the arguments, and if they are not set or have default values, we set a filter boolean to false.
        // Only when filtering needs to happen the filterrecord will be invoked, for performance.
        let mut filter: bool = false;
        let mut manipulate: bool = false;
        let _mmq = minmappingquality.unwrap_or(0);
        let _sfi = samflaginclude.unwrap_or(0);
        let _sfe = samflagexclude.unwrap_or(0);
        let _mifl = minfraglen.unwrap_or(0);
        let _mafl = maxfraglen.unwrap_or(0);
        let _mnase = mnase.unwrap_or(false);
        let _offset = offset.unwrap_or((0, 0));
        let _frs = filterrnastrand.unwrap_or(String::from("None"));
        let _extend = extendreads.unwrap_or(false);
        let _extendreadslen = extendreadslen.unwrap_or(0);
        let _extend_auto = _extend && _extendreadslen == 0;
        let _center = centerreads.unwrap_or(false);

        // Set the manipulate bool for a quick escape in case manipulation is not needed.
        if _offset != (0, 0) || _mnase || _extend || _center {
            manipulate = true;
        }

        // Set the filter bool for a quick escape in case filtering is not needed.
        if _mmq > 0
            || _sfi > 0
            || _sfe > 0
            || _mifl > 0
            || _mafl > 0
            || _frs != "None"
            || blacklist.is_some()
        {
            filter = true;
        }

        let blacklist_index = if blacklist.is_some() {
            Some(BlacklistIndex::from_regions(blacklist.as_ref().unwrap()))
        } else {
            None
        };

        Self {
            blacklist_index,
            blacklist: blacklist,
            minmappingquality: _mmq,
            samflaginclude: _sfi,
            samflagexclude: _sfe,
            minfraglen: _mifl,
            maxfraglen: _mafl,
            mnase: _mnase,
            offset: _offset,
            filterrnastrand: _frs,
            extendreads: _extend,
            extendreadslen: _extendreadslen,
            extendreads_auto: _extend_auto,
            centerreads: _center,
            filter: filter,
            manipulate: manipulate,
        }
    }

    pub fn set_extendreadslen(&mut self, bamfile: &str, nproc: usize, regions: &Vec<Region>) {
        const FREAD: u16 = 0x40;
        let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();
        let fraglens: Vec<u32> = pool.install(|| {
            regions
                .par_iter()
                .flat_map(|i| {
                    let mut bam = IndexedReader::from_path(bamfile).unwrap();
                    bam.fetch((i.chrom.as_str(), i.get_startu(), i.get_endu()))
                        .expect(&format!("Error fetching region: {:?}", i));
                    let mut fraglens: Vec<u32> = vec![];
                    for record in bam.records() {
                        let record = record.expect("Error parsing record.");
                        if record.is_paired()
                            && record.is_proper_pair()
                            && (record.flags() & FREAD != 0)
                        {
                            fraglens.push(record.insert_size().abs() as u32);
                        }
                    }
                    fraglens
                })
                .collect()
        });
        if fraglens.len() > 0 {
            let medlen = median(fraglens);
            self.extendreadslen = medlen.ceil() as u32;
        } else {
            panic!("No proper pairs found in the given regions. Please check your input.");
        }
    }

    pub fn filter_record(&self, rec: &Record, chrom: &str) -> bool {
        // Decides filtering of a record. The bool return is used to 'continue', i.e. skip the record.
        if rec.is_unmapped() {
            return true;
        } else if !self.filter {
            return false;
        } else {
            // True filtering.
            // quality > samflags > min/max fraglen
            // quality
            if rec.mapq() < self.minmappingquality {
                return true;
            }
            // samflags
            if self.samflaginclude > 0 {
                if (rec.flags() & self.samflaginclude) != self.samflaginclude {
                    return true;
                }
            }
            if self.samflagexclude > 0 {
                if (rec.flags() & self.samflagexclude) != 0 {
                    return true;
                }
            }
            // min/max fraglen
            if self.minfraglen != 0 || self.maxfraglen != 0 {
                if rec.is_paired() {
                    if rec.insert_size().abs() < self.minfraglen as i64
                        || rec.insert_size().abs() > self.maxfraglen as i64
                    {
                        return true;
                    }
                } else {
                    let fragsize: u32 = rec
                        .aligned_blocks()
                        .map(|x| x[1] as u32 - x[0] as u32)
                        .sum();
                    if fragsize < self.minfraglen || fragsize > self.maxfraglen {
                        return true;
                    }
                }
            }
            // filterrnastrand
            if self.filterrnastrand.as_str() != "None" {
                match (self.filterrnastrand.as_str(), rec.is_paired()) {
                    ("forward", true) => {
                        if !((rec.flags() & 144 == 128) || (rec.flags() & 96 == 64)) {
                            return true;
                        }
                    }
                    ("forward", false) => {
                        if !(rec.flags() & 16 == 16) {
                            return true;
                        }
                    }
                    ("reverse", true) => {
                        if !((rec.flags() & 144 == 144) || (rec.flags() & 96 == 96)) {
                            return true;
                        }
                    }
                    ("reverse", false) => {
                        if !(rec.flags() & 16 == 0) {
                            return true;
                        }
                    }
                    _ => {
                        panic!(
                            "filterrnastrand should be either forward or reverse. {:?} is not supported.",
                            self.filterrnastrand
                        )
                    }
                }
            }
            if self.blacklist.is_some() {
                return self.rec_in_blacklist(rec, chrom);
            }
            false
        }
    }

    pub fn manipulate_record(&self, rec: &Record) -> Option<Vec<u32>> {
        // Not just simple yes - no filtering, but actually manipulating the record.
        // In general, this is the case for MNase mode, offset, extendreads and centerreads.
        if self.mnase {
            // MNase mode, take only middle bps of the fragment

            // only retain proper pairs and forward read.
            let rinsertsize = rec.insert_size().abs() as u32;
            if rec.is_proper_pair() && !rec.is_reverse() && rinsertsize > 1 {
                // Not sure why the insert size test is here, but note we always filter prior to manipulating.
                let recpos: u32 = rec.pos() as u32;
                let frag_start = recpos - 1 + rinsertsize / 2;

                if rinsertsize % 2 == 0 {
                    return Some((frag_start..frag_start + 2).collect());
                } else {
                    return Some((frag_start..frag_start + 4).collect());
                }
            }
            return None;
        }
        if self.offset != (0, 0) {
            let mut blockvec: Vec<u32> = self.rec_blocks_for_offset(rec);

            let blocklen = blockvec.len() as i32;

            // Convert potential negative indices to positive indices
            // It could be that for the offset only one value is given, in which case we only use that site
            if self.offset.1 == 0 {
                let pos = if self.offset.0 < 0 {
                    blocklen + self.offset.0
                } else {
                    self.offset.0 - 1
                };
                if pos < 0 || pos >= blocklen {
                    return None;
                }

                if rec.is_reverse() {
                    blockvec.reverse();
                    blockvec = blockvec[pos as usize..pos as usize + 1].to_vec();
                    blockvec.reverse();
                    return Some(blockvec);
                } else {
                    blockvec = blockvec[pos as usize..pos as usize + 1].to_vec();
                    return Some(blockvec);
                }
            } else {
                let start = if self.offset.0 < 0 {
                    blocklen + self.offset.0
                } else {
                    self.offset.0 - 1
                };
                let end = if self.offset.1 < 0 {
                    blocklen + self.offset.1 + 1
                } else {
                    self.offset.1
                };

                if start < 0 || end < 0 || start >= blocklen || end > blocklen || start >= end {
                    return None;
                }
                if rec.is_reverse() {
                    blockvec.reverse();
                    blockvec = blockvec[start as usize..end as usize].to_vec();
                    blockvec.reverse();
                    return Some(blockvec);
                } else {
                    blockvec = blockvec[start as usize..end as usize].to_vec();
                    return Some(blockvec);
                }
            }
        }
        if self.extendreads {
            // Extend reads
            let blockvec = self.rec_extension(rec);
            return Some(blockvec);
        }
        return None;
    }

    fn is_proper_pair_precise(&self, rec: &Record, max_paired_fragment_length: i64) -> bool {
        if !rec.is_proper_pair() {
            return false;
        }
        if rec.tid() != rec.mtid() {
            return false;
        }
        if rec.insert_size().abs() > max_paired_fragment_length {
            return false;
        }
        if rec.is_reverse() == rec.is_mate_reverse() {
            return false;
        }
        if rec.is_reverse() {
            rec.pos() >= rec.mpos()
        } else {
            rec.pos() <= rec.mpos()
        }
    }

    pub fn rec_blocks_for_offset(&self, rec: &Record) -> Vec<u32> {
        let mut blocks: Vec<(i64, i64)> = rec
            .aligned_blocks()
            .map(|x| (x[0] as i64, x[1] as i64))
            .collect();

        if self.extendreads {
            let query_len: i64 = rec.seq_len_from_cigar(false) as i64;
            let max_paired_fragment_length: i64 = if self.maxfraglen > 0 {
                self.maxfraglen as i64
            } else {
                4 * self.extendreadslen as i64
            };
            let max_paired_fragment_length = if max_paired_fragment_length <= 0 {
                1000
            } else {
                max_paired_fragment_length
            };

            if self.is_proper_pair_precise(rec, max_paired_fragment_length) {
                if rec.is_reverse() {
                    let foo = (rec.mpos(), rec.reference_start());
                    if foo.0 < foo.1 {
                        blocks.insert(0, foo);
                    }
                } else {
                    let foo = (
                        rec.reference_end(),
                        rec.reference_end() + rec.insert_size().abs() - query_len,
                    );
                    if foo.0 < foo.1 {
                        blocks.push(foo);
                    }
                }
            } else {
                if rec.is_reverse() {
                    let mut start = rec.reference_start() - self.extendreadslen as i64 + query_len;
                    if start < 0 {
                        start = 0;
                    }
                    let foo = (start, rec.reference_start());
                    if foo.0 < foo.1 {
                        blocks.insert(0, foo);
                    }
                } else {
                    let foo = (
                        rec.reference_end(),
                        rec.reference_end() + self.extendreadslen as i64 - query_len,
                    );
                    if foo.0 < foo.1 {
                        blocks.push(foo);
                    }
                }
            }
        }

        let mut stretch: Vec<u32> = Vec::new();
        for (s, e) in blocks {
            stretch.extend((s.max(0) as u32)..(e.max(0) as u32));
        }

        stretch
    }

    pub fn rec_extension(&self, rec: &Record) -> Vec<u32> {
        let read_len: i64 = rec.seq_len_from_cigar(false) as i64;
        let is_pe = rec.is_proper_pair();

        if !self.extendreads_auto && (self.extendreadslen as i64) <= read_len {
            let mut blockvec: Vec<u32> = Vec::new();
            rec.aligned_blocks().for_each(|x| {
                blockvec.extend(x[0] as u32..x[1] as u32);
            });
            return blockvec;
        }

        let (mut fragment_start, mut fragment_end): (i64, i64) = if is_pe {
            if rec.is_reverse() {
                (rec.mpos(), rec.reference_end())
            } else {
                (
                    rec.reference_start(),
                    rec.reference_start() + rec.insert_size().abs(),
                )
            }
        } else if rec.is_reverse() {
            (
                rec.reference_end() - self.extendreadslen as i64,
                rec.reference_end(),
            )
        } else {
            (
                rec.reference_start(),
                rec.reference_start() + self.extendreadslen as i64,
            )
        };

        if self.centerreads && fragment_end > fragment_start {
            let fragment_center =
                fragment_end as f64 - (fragment_end - fragment_start) as f64 / 2.0;
            fragment_start = (fragment_center - read_len as f64 / 2.0).trunc() as i64;
            fragment_end = fragment_start + read_len;
        }

        (fragment_start.max(0) as u32..fragment_end.max(0) as u32).collect()
    }

    pub fn rec_in_blacklist(&self, rec: &Record, chrom: &str) -> bool {
        let pos = rec.pos() as u32;
        let idx = self.blacklist_index.as_ref().unwrap();
        if idx.contains(chrom, pos) {
            return true;
        }
        let end = rec.seq_len_from_cigar(false) as u32 + pos;
        if idx.contains(chrom, end) {
            return true;
        }
        false
    }
}
