use crate::normalization::*;

mod scale_factor_tests {
    use super::*;

    fn sf(method: Normalization, preset: Option<f32>) -> Result<f32, NormalizationError> {
        scale_factor(method, 1000, 10, 10000, 20.0, preset)
    }

    #[test]
    fn test_scale_factor_rpkm() {
        assert_eq!(sf(Normalization::Rpkm, None).unwrap(), 99999.99);
    }

    #[test]
    fn test_scale_factor_cpm() {
        assert_eq!(sf(Normalization::Cpm, None).unwrap(), 999.99994);
    }

    #[test]
    fn test_scale_factor_bpm() {
        assert_eq!(sf(Normalization::Bpm, None).unwrap(), 999.99994);
    }

    #[test]
    fn test_scale_factor_rpgc() {
        assert_eq!(sf(Normalization::Rpgc, None).unwrap(), 0.5);
    }

    #[test]
    fn test_scale_factor_none() {
        assert_eq!(sf(Normalization::None, None).unwrap(), 1.0);
    }

    #[test]
    fn test_preset_scalefactor() {
        assert_eq!(sf(Normalization::Rpkm, Some(5.5)).unwrap(), 5.5);
    }

    #[test]
    fn test_rpgc_needs_genome_size() {
        assert_eq!(
            scale_factor(Normalization::Rpgc, 1000, 10, 0, 20.0, None),
            Err(NormalizationError::MissingEffectiveGenomeSize)
        );
    }

    #[test]
    fn test_preset_skips_rpgc_genome_size_check() {
        assert_eq!(
            scale_factor(Normalization::Rpgc, 1000, 10, 0, 20.0, Some(2.0)),
            Ok(2.0)
        );
    }
}

mod parse_tests {
    use super::*;

    #[test]
    fn test_parse_normalization() {
        assert_eq!("RPGC".parse::<Normalization>(), Ok(Normalization::Rpgc));
        assert_eq!("None".parse::<Normalization>(), Ok(Normalization::None));
        assert!(matches!(
            "rpkm".parse::<Normalization>(),
            Err(NormalizationError::UnknownMethod { .. })
        ));
    }

    #[test]
    fn test_parse_scale_factors_method() {
        assert_eq!("readCount".parse::<ScaleFactorsMethod>(), Ok(ScaleFactorsMethod::ReadCount));
        assert!(matches!(
            "bogus".parse::<ScaleFactorsMethod>(),
            Err(NormalizationError::UnknownMethod { .. })
        ));
    }
}

mod bamcompare_scale_factor_tests {
    use super::*;

    fn samples() -> [SampleStats; 2] {
        [
            SampleStats { mapped: 1000, readlen: 20.0 },
            SampleStats { mapped: 2000, readlen: 20.0 },
        ]
    }

    #[test]
    fn test_bamcompare_readcount() {
        let (sf1, sf2) = scale_factor_bamcompare(
            ScaleFactorsMethod::ReadCount,
            samples(),
            10,
            10000,
            Normalization::Rpkm,
        )
        .unwrap();
        assert_eq!(sf1, 1.0);
        assert_eq!(sf2, 0.5);
    }

    #[test]
    fn test_bamcompare_independent_norm() {
        let (sf1, sf2) = scale_factor_bamcompare(
            ScaleFactorsMethod::None,
            samples(),
            10,
            10000,
            Normalization::Cpm,
        )
        .unwrap();
        assert_eq!(sf1, 999.99994);
        assert!((sf2 - 500.0).abs() < 1e-3, "sf2 = {sf2}");
    }

    #[test]
    fn test_bamcompare_rpgc_needs_genome_size() {
        assert_eq!(
            scale_factor_bamcompare(
                ScaleFactorsMethod::None,
                samples(),
                10,
                0,
                Normalization::Rpgc,
            ),
            Err(NormalizationError::MissingEffectiveGenomeSize)
        );
    }
}
