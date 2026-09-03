use crate::covcalc::{Region, Revalue};
use crate::filtering::{Alignmentfilters, BlacklistIndex};

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
