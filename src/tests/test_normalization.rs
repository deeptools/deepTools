use crate::normalization::*;

mod scale_factor_tests {
    use super::*;

    #[test]
    fn test_scale_factor_rpkm() {
        let nm = "RPKM";
        let sf = scale_factor(
            nm,
            1000, 10, 10000,
            20.0, 20.0,
            1.0, &false
        );
        assert_eq!(sf, 99999.99);
    }

    #[test]
    fn test_scale_factor_cpm() {
        let nm = "CPM";
        let sf = scale_factor(
            nm,
            1000, 10, 10000,
            20.0, 20.0,
            1.0, &false
        );
        assert_eq!(sf, 999.99994);
    }

    #[test]
    fn test_scale_factor_bpm() {
        let nm = "BPM";
        let sf = scale_factor(
            nm,
            1000, 10, 10000,
            20.0, 20.0,
            1.0, &false
        );
        assert_eq!(sf, 999.99994);
    }

    #[test]
    fn test_scale_factor_rpgc() {
        let nm = "RPGC";
        let sf = scale_factor(
            nm,
            1000, 10, 10000,
            20.0, 20.0,
            1.0, &false
        );
        assert_eq!(sf, 0.5);
    }

    #[test]
    fn test_preset_scalefactor() {
        let nm = "RPKM";
        let sf = scale_factor(
            nm,
            1000, 10, 10000,
            20.0, 20.0,
            5.5, &false
        );
        assert_eq!(sf, 5.5);
    }
}

mod bamcompare_scale_factor_tests{
    use super::*;
    
    #[test]
    fn test_bamcompare_readcount() {
        let nm = "readCount";
        let (sf1, sf2) = scale_factor_bamcompare(
            nm,
            1000, 2000,
            10, 10000,
            "RPKM", 
            20.0, 20.0,
            1.0, 1.0
        );
        assert_eq!(sf1, 1.0);
        assert_eq!(sf2, 0.5);
    }
}