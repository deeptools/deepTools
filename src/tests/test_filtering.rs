use crate::alignmentsieve::with_position;
use crate::covcalc::{Region, Revalue};
use crate::filtering::{Alignmentfilters, BlacklistIndex, center_window, fragment_window};
use rust_htslib::bam::{
    self,
    record::{Aux, Cigar, CigarString},
};

fn region(chrom: &str, start: u32, end: u32) -> Region {
    Region {
        chrom: chrom.to_string(),
        start: Revalue::U(start),
        end: Revalue::U(end),
        score: String::from("."),
        strand: String::from("+"),
        name: String::from("test"),
        regionlength: end - start,
    }
}

mod blacklist_index_tests {
    use super::*;

    #[test]
    fn test_contains_inside_region() {
        let regions = vec![region("chr1", 100, 200)];
        let idx = BlacklistIndex::from_regions(&regions);
        assert!(idx.contains("chr1", 150));
    }

    #[test]
    fn test_contains_at_region_boundaries() {
        let regions = vec![region("chr1", 100, 200)];
        let idx = BlacklistIndex::from_regions(&regions);
        assert!(idx.contains("chr1", 100));
        assert!(idx.contains("chr1", 200));
    }

    #[test]
    fn test_does_not_contain_outside_region() {
        let regions = vec![region("chr1", 100, 200)];
        let idx = BlacklistIndex::from_regions(&regions);
        assert!(!idx.contains("chr1", 99));
        assert!(!idx.contains("chr1", 201));
    }

    #[test]
    fn test_unknown_chrom_returns_false() {
        let regions = vec![region("chr1", 100, 200)];
        let idx = BlacklistIndex::from_regions(&regions);
        assert!(!idx.contains("chr2", 150));
    }

    #[test]
    fn test_multiple_non_overlapping_regions() {
        let regions = vec![region("chr1", 100, 200), region("chr1", 500, 600)];
        let idx = BlacklistIndex::from_regions(&regions);
        assert!(idx.contains("chr1", 550));
        assert!(!idx.contains("chr1", 300));
    }

    #[test]
    fn test_unsorted_input_regions_are_handled() {
        let regions = vec![region("chr1", 500, 600), region("chr1", 100, 200)];
        let idx = BlacklistIndex::from_regions(&regions);
        assert!(idx.contains("chr1", 150));
        assert!(idx.contains("chr1", 550));
    }

    #[test]
    fn test_empty_index_never_contains() {
        let idx = BlacklistIndex::from_regions(&[]);
        assert!(!idx.contains("chr1", 1));
    }
}

mod alignmentfilters_new_tests {
    use super::*;

    #[test]
    fn test_defaults_disable_filter_and_manipulate() {
        let af = Alignmentfilters::new(
            None, None, None, None, None, None, None, None, None, None, None, None,
        );
        assert!(!af.filter);
        assert!(!af.manipulate);
    }

    #[test]
    fn test_min_mapping_quality_enables_filter() {
        let af = Alignmentfilters::new(
            None,
            Some(10),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        );
        assert!(af.filter);
        assert!(!af.manipulate);
        assert_eq!(af.minmappingquality, 10);
    }

    #[test]
    fn test_blacklist_enables_filter() {
        let regions = vec![region("chr1", 100, 200)];
        let af = Alignmentfilters::new(
            Some(regions),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        );
        assert!(af.filter);
    }

    #[test]
    fn test_extendreads_enables_manipulate_and_defaults_to_auto() {
        let af = Alignmentfilters::new(
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            Some(true),
            None,
            None,
        );
        assert!(af.manipulate);
        assert!(af.extendreads);
        assert!(af.extendreads_auto);
    }

    #[test]
    fn test_extendreads_with_explicit_len_disables_auto() {
        let af = Alignmentfilters::new(
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            Some(true),
            Some(150),
            None,
        );
        assert!(af.manipulate);
        assert!(!af.extendreads_auto);
        assert_eq!(af.extendreadslen, 150);
    }

    #[test]
    fn test_offset_enables_manipulate_only() {
        let af = Alignmentfilters::new(
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            Some((1, 2)),
            None,
            None,
            None,
            None,
        );
        assert!(af.manipulate);
        assert!(!af.filter);
    }

    #[test]
    fn test_mnase_enables_manipulate_only() {
        let af = Alignmentfilters::new(
            None,
            None,
            None,
            None,
            None,
            None,
            Some(true),
            None,
            None,
            None,
            None,
            None,
        );
        assert!(af.manipulate);
        assert!(!af.filter);
    }

    #[test]
    fn test_centerreads_alone_enables_manipulate() {
        let af = Alignmentfilters::new(
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            Some(true),
        );
        assert!(af.manipulate);
    }

    #[test]
    fn test_filterrnastrand_enables_filter() {
        let af = Alignmentfilters::new(
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            Some(String::from("forward")),
            None,
            None,
            None,
        );
        assert!(af.filter);
    }
}

mod mnase_centre_tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString, Record};

    fn mnase_filters() -> Alignmentfilters {
        Alignmentfilters::new(
            None, None, None, None, None, None,
            Some(true), // MNase mode
            None, None, None, None, None,
        )
    }

    fn proper_pair_forward_read(pos: i64, insert_size: i64) -> Record {
        let mut rec = Record::new();
        rec.set(
            b"frag",
            Some(&CigarString(vec![Cigar::Match(50)])),
            &[b'A'; 50],
            &[30u8; 50],
        );
        rec.set_tid(0);
        rec.set_pos(pos);
        rec.set_flags(99); // paired, proper pair, mate reverse, first in pair
        rec.set_mtid(0);
        rec.set_mpos(pos + insert_size - 50);
        rec.set_insert_size(insert_size);
        rec
    }

    #[test]
    fn test_mnase_even_fragment_two_central_bases() {
        // fragment [100, 250): central bases 174 and 175
        let af = mnase_filters();
        let rec = proper_pair_forward_read(100, 150);
        assert_eq!(af.manipulate_record(&rec), Some(vec![174, 175]));
    }

    #[test]
    fn test_mnase_odd_fragment_three_central_bases() {
        // fragment [100, 249): centre 174, documented bases 173, 174, 175
        let af = mnase_filters();
        let rec = proper_pair_forward_read(100, 149);
        assert_eq!(af.manipulate_record(&rec), Some(vec![173, 174, 175]));
        // fragment [1000, 1131): centre 1065
        let rec = proper_pair_forward_read(1000, 131);
        assert_eq!(af.manipulate_record(&rec), Some(vec![1064, 1065, 1066]));
    }
}

fn bam_rec(paired: bool, proper: bool, rev: bool, pos: i64, tlen: i64, mpos: i64) -> bam::Record {
    let mut r = bam::Record::new();
    let cigar = CigarString(vec![Cigar::Match(50)]);
    r.set(b"q", Some(&cigar), &vec![b'A'; 50], &vec![30u8; 50]);
    r.set_pos(pos);
    r.set_tid(0);
    let mut flags = 0u16;
    if paired {
        flags |= 0x1;
    }
    if proper {
        flags |= 0x2;
    }
    if rev {
        flags |= 0x10;
    }
    r.set_flags(flags);
    r.set_insert_size(tlen);
    r.set_mpos(mpos);
    r
}

#[test]
fn test_fragment_window_single_end() {
    let fwd = bam_rec(false, false, false, 1000, 0, -1);
    assert_eq!(fragment_window(&fwd, 200), (1000, 1200));
    let rev = bam_rec(false, false, true, 1000, 0, -1);
    assert_eq!(fragment_window(&rev, 200), (850, 1050));
    assert_eq!(fragment_window(&rev, 50), (1000, 1050));
}

#[test]
fn test_fragment_window_proper_pair_ignores_ext_len() {
    let fwd = bam_rec(true, true, false, 1000, 300, 1250);
    assert_eq!(fragment_window(&fwd, 50), (1000, 1300));
    // Reverse mate: tlen is negative, window runs from the mate start to this read's end.
    let rev = bam_rec(true, true, true, 1250, -300, 1000);
    assert_eq!(fragment_window(&rev, 50), (1000, 1300));
}

#[test]
fn test_fragment_window_paired_but_not_proper_is_treated_as_single_end() {
    // Same rule as Alignmentfilters::rec_extension (is_proper_pair).
    let fwd = bam_rec(true, false, false, 1000, 300, 1250);
    assert_eq!(fragment_window(&fwd, 200), (1000, 1200));
}

#[test]
fn test_center_window_centres_a_read_in_the_fragment() {
    assert_eq!(center_window(1000, 1300, 50), (1125, 1175));
    // Odd lengths truncate toward zero, matching rec_extension.
    assert_eq!(center_window(1000, 1301, 50), (1125, 1175));
}

#[test]
fn test_center_window_leaves_empty_window_alone() {
    assert_eq!(center_window(1000, 1000, 50), (1000, 1000));
}

#[test]
fn test_with_position_moves_read_keeps_seq_qual_and_tags() {
    let mut r = bam_rec(false, false, false, 1000, 0, -1);
    r.push_aux(b"CB", Aux::String("AAAC-1")).unwrap();
    let moved = with_position(&r, 2000, 50);
    assert_eq!(moved.pos(), 2000);
    assert_eq!(moved.cigar().end_pos(), 2050);
    assert_eq!(moved.seq().as_bytes(), r.seq().as_bytes());
    assert_eq!(moved.qual(), r.qual());
    match moved.aux(b"CB") {
        Ok(Aux::String(s)) => assert_eq!(s, "AAAC-1"),
        other => panic!("CB tag lost or changed: {other:?}"),
    }
}

#[test]
fn test_with_position_replaces_complex_cigar_with_single_match() {
    let mut r = bam::Record::new();
    let cigar = CigarString(vec![Cigar::Match(20), Cigar::Del(5), Cigar::Match(30)]);
    r.set(b"q", Some(&cigar), &vec![b'A'; 50], &vec![30u8; 50]);
    r.set_pos(100);
    let moved = with_position(&r, 500, 50);
    assert_eq!(moved.cigar().end_pos(), 550);
    assert_eq!(moved.cigar().len(), 1);
}
