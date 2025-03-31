use rust_htslib::bam::{Record, ext::BamRecordExtensions};
use crate::covcalc::Region;

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
    pub extendreads: u32,
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
        extendreads: Option<u32>,
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
        let _extend = extendreads.unwrap_or(0);
        let _center = centerreads.unwrap_or(false);

        // Set the manipulate bool for a quick escape in case manipulation is not needed.
        if _offset != (0, 0) || _mnase || _extend > 0 || _center {
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
            centerreads: _center,
            filter: filter,
            manipulate: manipulate,
        }
    }
    pub fn filter_record(&self, rec: &Record) -> bool {
        // Decides filtering of a record. The bool return is used to 'continue', i.e. skip the record.
        if rec.is_unmapped() {
            return true;
        } else if !self.filter {
            return false;
        } else {
            // True filtering.
            // quality > samflags > min/max fraglen
            // quality
            let mut skip: bool = false;
            if rec.mapq() < self.minmappingquality {
                skip = true;
            }
            // samflags
            if self.samflaginclude > 0 {
                if (rec.flags() & self.samflaginclude) == 0 {
                    skip = true;
                }
            }
            if self.samflagexclude > 0 {
                if (rec.flags() & self.samflagexclude) != 0 {
                    skip = true;
                }
            }
            // min/max fraglen
            if self.minfraglen != 0 || self.maxfraglen != 0 {
                if rec.is_paired() {
                    if rec.insert_size().abs() < self.minfraglen as i64 || rec.insert_size().abs() > self.maxfraglen as i64 {
                        skip = true;
                    }
                } else {
                    let fragsize: u32 = rec
                        .aligned_blocks()
                        .map(|x| x[1] as u32 - x[0] as u32)
                        .sum();
                    if fragsize < self.minfraglen || fragsize > self.maxfraglen {
                        skip = true;
                    }
                }
            }
            // filterrnastrand
            if self.filterrnastrand.as_str() != "None" {
                match (self.filterrnastrand.as_str(), rec.is_paired()) {
                    ("forward", true) => {
                        if !((rec.flags() & 144 == 128) || (rec.flags() & 96 == 64)) {
                            skip = true;
                        }
                    },
                    ("forward", false) => {
                        if !(rec.flags() & 16 == 16) {
                            skip = true;
                        }
                    },
                    ("reverse", true) => {
                        if !((rec.flags() & 144 == 144) || (rec.flags() & 96 == 96)) {
                            skip = true;
                        }
                    },
                    ("reverse", false) => {
                        if !(rec.flags() & 16 == 0) {
                            skip = true;
                        }
                    },
                    _ => {
                        panic!("filterrnastrand should be either forward or reverse. {:?} is not supported.", self.filterrnastrand)
                    },
                }
            }
            return skip;
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
            // Collect blocks and flatten them out.
            let mut blockvec: Vec<u32> = rec
                .aligned_blocks()
                .flat_map(|x| x[0] as u32..x[1] as u32)
                .collect();
            let blocklen = blockvec.len() as i32;

            // Convert potential negative indices to positive indices
            // It could be that for the offset only one value is given, in which case we only use that site
            if self.offset.1 == 0 {
                // 
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
        if self.extendreads > 0 {
            println!("extendreads implementation");
        }
        return None;
    }
}