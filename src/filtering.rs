use rust_htslib::bam::{Record, ext::BamRecordExtensions, IndexedReader, Read};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use crate::covcalc::Region;
use crate::calc::median;

#[derive(Clone)]
pub struct Alignmentfilters {
    pub blacklist: Option<Vec<Region>>,
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
        centerreads: Option<bool>
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
        let _offset =  offset.unwrap_or((0, 0));
        let _frs = filterrnastrand.unwrap_or(String::from("None"));
        let _extend = extendreads.unwrap_or(false);
        let _extendreadslen = extendreadslen.unwrap_or(0);
        let _center = centerreads.unwrap_or(false);

        // Set the manipulate bool for a quick escape in case manipulation is not needed.
        if _offset != (0, 0) || _mnase || _extend || _center {
            manipulate = true;
        }

        // Set the filter bool for a quick escape in case filtering is not needed.
        if _mmq > 0 || _sfi > 0 || _sfe > 0 || _mifl > 0 || _mafl > 0 || _frs != "None" || blacklist.is_some() {
            filter = true;
        }

        Self {
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
            centerreads: _center,
            filter: filter,
            manipulate: manipulate,
        }
    }

    pub fn set_extendreadslen(&mut self, bamfile: &str, nproc: usize, regions: &Vec<Region> ) {
        const FREAD: u16 = 0x40;
        let pool = ThreadPoolBuilder::new().num_threads(nproc).build().unwrap();
        let fraglens: Vec<u32> = pool.install(|| {
            regions
                .par_iter()
                .flat_map(|i| {
                    let mut bam = IndexedReader::from_path(bamfile).unwrap();
                    bam.fetch((i.chrom.as_str(), i.get_startu(), i.get_endu()) )
                        .expect(&format!("Error fetching region: {:?}", i));
                    let mut fraglens: Vec<u32> = vec![];
                    for record in bam.records() {
                        let record = record.expect("Error parsing record.");
                        if record.is_paired() && record.is_proper_pair() && (record.flags() & FREAD != 0) {
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
                    if rec.insert_size().abs() < self.minfraglen as i64 || rec.insert_size().abs() > self.maxfraglen as i64 {
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
                    },
                    ("forward", false) => {
                        if !(rec.flags() & 16 == 16) {
                            return true;
                        }
                    },
                    ("reverse", true) => {
                        if !((rec.flags() & 144 == 144) || (rec.flags() & 96 == 96)) {
                            return true;
                        }
                    },
                    ("reverse", false) => {
                        if !(rec.flags() & 16 == 0) {
                            return true;
                        }
                    },
                    _ => {
                        panic!("filterrnastrand should be either forward or reverse. {:?} is not supported.", self.filterrnastrand)
                    },
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
                    return Some(
                        (frag_start..frag_start + 2).collect()
                    );
                } else {
                    return Some(
                        (frag_start..frag_start+4).collect()
                    );
                }
            }
            return None;
        }
        if self.offset != (0, 0) {

            let mut blockvec: Vec<u32> = if self.extendreads {
                self.rec_extension(rec)
            } else {
                rec
                    .aligned_blocks()
                    .flat_map(|x| x[0] as u32..x[1] as u32)
                    .collect()
            };

            let blocklen = blockvec.len() as i32;

            // Convert potential negative indices to positive indices
            // It could be that for the offset only one value is given, in which case we only use that site
            if self.offset.1 == 0 {
                let pos = if self.offset.0 < 0 {blocklen + self.offset.0 } else {self.offset.0 - 1};
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
                let start = if self.offset.0 < 0 { blocklen + self.offset.0 } else { self.offset.0 -1};
                let end = if self.offset.1 < 0 { blocklen + self.offset.1 + 1 } else { self.offset.1 };

                // if the range falls outside the vec, return none (retain deeptools 3 behavior)
                if start < 0 || end < 0 || start >= blocklen || end >= blocklen || start >= end {
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
            return Some(blockvec)
        }
        return None;
    }

    pub fn rec_extension(&self, rec: &Record) -> Vec<u32> {
        // extend the reads. if the read is a proper pair, we get the fragment length from there.
        // If not (or if se), then extendreadslen is used.
        // Note that extendsreadslen is already populated at this stage, either by CLI (se) or calculated (pe).
        let mut blockvec: Vec<u32> = Vec::new();
        let mut blocklen: u32 = 0;

        rec
            .aligned_blocks()
            .for_each(|x| {
                let _s = x[0] as u32;
                let _e = x[1] as u32;
                blockvec.extend(_s.._e);
                blocklen += _e - _s;
            });

        if rec.is_proper_pair() {
            // Proper pairs
            if rec.is_reverse() {
                let ns = rec.mpos() as u32;
                let ne = rec.reference_start() as u32;
                if ns < ne {
                    // blockvec.splice(0..0,ns..ne);
                    let mut new_blockvec = Vec::with_capacity((ne - ns) as usize + blockvec.len());
                    new_blockvec.extend(ns..ne);
                    new_blockvec.extend(blockvec);
                    blockvec = new_blockvec;
                }
            } else {
                let ns = rec.reference_end() as u32;
                let ne: u32 = ns + rec.insert_size().abs() as u32 - rec.seq_len_from_cigar(false) as u32;
                if ns < ne {
                    blockvec.extend(ns..ne);
                }
            }
        } else {
            // non proper pairs -> 'se mode'
            if rec.is_reverse() {
                let ns: u32;
                let _rem = self.extendreadslen - rec.seq_len_from_cigar(false) as u32;
                if _rem > rec.reference_start() as u32 {
                    ns = 0;
                } else {
                    ns = rec.reference_start() as u32 - _rem;
                }
                let ne = rec.reference_start() as u32;
                if ns < ne {
                    //blockvec.splice(0..0,ns..ne);
                    let mut new_blockvec = Vec::with_capacity((ne - ns) as usize + blockvec.len());
                    new_blockvec.extend(ns..ne);
                    new_blockvec.extend(blockvec);
                    blockvec = new_blockvec;
                }
            } else {
                let ns = rec.reference_end() as u32;
                let ne: u32 = ns + self.extendreadslen - rec.seq_len_from_cigar(false) as u32;
                if ns < ne {
                    blockvec.extend(ns..ne );
                }
            }
        }
        if self.centerreads {
            let centerpoint = (blockvec.len() as u32 - blocklen) / 2;
            return blockvec[centerpoint as usize..(centerpoint + blocklen) as usize].to_vec();
        }
        return blockvec;
    }

    pub fn rec_in_blacklist(&self, rec: &Record, chrom: &str ) -> bool {
        for region in self.blacklist.as_ref().unwrap().iter() {
            if region.chrom == chrom {
                let pos = rec.pos() as u32;
                if region.get_startu() <= pos as u32 && pos as u32 <= region.get_endu() {
                    return true;
                }
                let end = rec.seq_len_from_cigar(false) as u32 + pos;
                if region.get_startu() <= end && end <= region.get_endu() {
                    return true;
                }
            }
        }
        false
    }

}