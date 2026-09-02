use crate::covcalc::{Region, Revalue};

fn region_u(strand: &str, start: u32, end: u32) -> Region {
    Region {
        chrom: String::from("chr1"),
        start: Revalue::U(start),
        end: Revalue::U(end),
        score: String::from("."),
        strand: strand.to_string(),
        name: String::from("test"),
        regionlength: end - start,
    }
}

fn region_v(strand: &str, starts: Vec<u32>, ends: Vec<u32>) -> Region {
    let regionlength = starts.iter().zip(ends.iter()).map(|(s, e)| e - s).sum();
    Region {
        chrom: String::from("chr1"),
        start: Revalue::V(starts),
        end: Revalue::V(ends),
        score: String::from("."),
        strand: strand.to_string(),
        name: String::from("test"),
        regionlength,
    }
}

mod get_startu_endu_tests {
    use super::*;

    #[test]
    fn test_startu_endu_non_metagene() {
        let r = region_u("+", 100, 200);
        assert_eq!(r.get_startu(), 100);
        assert_eq!(r.get_endu(), 200);
    }

    #[test]
    fn test_startu_endu_metagene_uses_outer_bounds() {
        let r = region_v("+", vec![100, 150], vec![120, 200]);
        assert_eq!(r.get_startu(), 100);
        assert_eq!(r.get_endu(), 200);
    }
}

mod anchorpoint_tss_tests {
    use super::*;

    #[test]
    fn test_tss_plus_strand_non_metagene() {
        let r = region_u("+", 100, 200);
        assert_eq!(r.get_anchorpoint("TSS"), 100);
    }

    #[test]
    fn test_tss_dot_strand_treated_as_plus() {
        let r = region_u(".", 100, 200);
        assert_eq!(r.get_anchorpoint("TSS"), 100);
    }

    #[test]
    fn test_tss_minus_strand_non_metagene() {
        let r = region_u("-", 100, 200);
        assert_eq!(r.get_anchorpoint("TSS"), 200);
    }

    #[test]
    fn test_tss_plus_strand_metagene_uses_first_exon_start() {
        let r = region_v("+", vec![100, 150], vec![120, 200]);
        assert_eq!(r.get_anchorpoint("TSS"), 100);
    }

    #[test]
    fn test_tss_minus_strand_metagene_uses_last_exon_end() {
        let r = region_v("-", vec![100, 150], vec![120, 200]);
        assert_eq!(r.get_anchorpoint("TSS"), 200);
    }

    #[test]
    #[should_panic(expected = "Strand should either be")]
    fn test_tss_invalid_strand_panics() {
        let r = region_u("x", 100, 200);
        r.get_anchorpoint("TSS");
    }
}

mod anchorpoint_tes_tests {
    use super::*;

    #[test]
    fn test_tes_plus_strand_non_metagene() {
        let r = region_u("+", 100, 200);
        assert_eq!(r.get_anchorpoint("TES"), 200);
    }

    #[test]
    fn test_tes_minus_strand_non_metagene() {
        let r = region_u("-", 100, 200);
        assert_eq!(r.get_anchorpoint("TES"), 100);
    }

    #[test]
    fn test_tes_plus_strand_metagene_uses_last_exon_end() {
        let r = region_v("+", vec![100, 150], vec![120, 200]);
        assert_eq!(r.get_anchorpoint("TES"), 200);
    }

    #[test]
    fn test_tes_minus_strand_metagene_uses_first_exon_start() {
        let r = region_v("-", vec![100, 150], vec![120, 200]);
        assert_eq!(r.get_anchorpoint("TES"), 100);
    }
}

mod anchorpoint_center_tests {
    use super::*;

    #[test]
    fn test_center_non_metagene_is_midpoint() {
        let r = region_u("+", 100, 200);
        assert_eq!(r.get_anchorpoint("center"), 150);
    }

    #[test]
    fn test_center_ignores_strand() {
        let plus = region_u("+", 100, 200);
        let minus = region_u("-", 100, 200);
        assert_eq!(
            plus.get_anchorpoint("center"),
            minus.get_anchorpoint("center")
        );
    }

    #[test]
    fn test_center_metagene_walks_exon_blocks() {
        let r = region_v("+", vec![100, 150], vec![120, 200]);
        assert_eq!(r.get_anchorpoint("center"), 165);
    }

    #[test]
    fn test_center_metagene_midpoint_in_first_exon() {
        let r = region_v("+", vec![100, 300], vec![180, 320]);
        assert_eq!(r.get_anchorpoint("center"), 150);
    }
}

mod anchorpoint_invalid_referencepoint_tests {
    use super::*;

    #[test]
    #[should_panic(expected = "Reference should either be TSS, TES or center")]
    fn test_unknown_referencepoint_panics() {
        let r = region_u("+", 100, 200);
        r.get_anchorpoint("bogus");
    }
}
