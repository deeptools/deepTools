use crate::calc::*;
use ndarray::{Array1, Array2};

mod vector_calculations_tests {
    use super::*;

    #[test]
    fn test_median() {
        let v: Vec<u32> = vec![1, 2, 3, 4, 5];
        assert_eq!(median(v), 3.0);
    }

    #[test]
    fn test_mean_float() {
        let v: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(mean_float(&v.iter().collect::<Vec<_>>()), 3.0);
    }

    #[test]
    fn test_median_float() {
        let v: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(median_float(&v.iter().collect::<Vec<_>>()), 3.0);
    }

    #[test]
    fn test_min_float() {
        let v: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(min_float(&v.iter().collect::<Vec<_>>()), 1.0);
    }

    #[test]
    fn test_max_float() {
        let v: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(max_float(&v.iter().collect::<Vec<_>>()), 5.0);
    }

    #[test]
    fn test_sum_float() {
        let v: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(sum_float(&v.iter().collect::<Vec<_>>()), 15.0);
    }

    #[test]
    fn test_std_float() {
        let v: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(std_float(&v.iter().collect::<Vec<_>>()), 1.4142135);

        let v: Vec<f32> = vec![1.0, 1.0, 1.0];
        assert_eq!(std_float(&v.iter().collect::<Vec<_>>()), 0.0);

        assert_eq!(std_float(&Vec::<&f32>::new()), 0.0);
    }
}

mod scalefactor_calculations_tests {
    use super::*;
    #[test]
    fn test_calc_ratio() {
        let cov1: f32 = 5.0;
        let cov2: f32 = 10.0;

        let sf1: f32 = 1.0;
        let sf2: f32 = 2.0;

        let pc1: f32 = 1.0;
        let pc2: f32 = 2.0;

        let log2 = "log2";
        let ratio = "ratio";
        let reciprocal_ratio = "reciprocal_ratio";
        let subtract = "subtract";

        let r1 = calc_ratio(cov1, cov2, &sf1, &sf2, &pc1, &pc2, log2);
        let r2 = calc_ratio(cov1, cov2, &sf1, &sf2, &pc1, &pc2, ratio);
        let r3 = calc_ratio(cov1, cov2, &sf1, &sf2, &pc1, &pc2, reciprocal_ratio);
        let r4 = calc_ratio(cov1, cov2, &sf1, &sf2, &pc1, &pc2, subtract);
        assert_eq!(r1, -1.87);
        assert_eq!(r2, 0.27);
        // 6/22 < 1, so the reciprocal ratio is -22/6
        assert_eq!(r3, -3.67);
        assert_eq!(r4, -16.0);
    }

    #[test]
    fn test_calc_ratio_reciprocal_ratio_sign_and_orientation() {
        // a/b if a/b >= 1 else -b/a (the 3.5.x getRatio doctest values)
        let one: f32 = 1.0;
        let zero: f32 = 0.0;
        assert_eq!(calc_ratio(2.0, 1.0, &one, &one, &zero, &zero, "reciprocal_ratio"), 2.0);
        assert_eq!(calc_ratio(1.0, 2.0, &one, &one, &zero, &zero, "reciprocal_ratio"), -2.0);
        assert_eq!(calc_ratio(1.0, 1.0, &one, &one, &zero, &zero, "reciprocal_ratio"), 1.0);
        assert_eq!(calc_ratio(3.0, 2.0, &one, &one, &zero, &zero, "reciprocal_ratio"), 1.5);
        assert_eq!(calc_ratio(2.0, 3.0, &one, &one, &zero, &zero, "reciprocal_ratio"), -1.5);
    }

    #[test]
    fn test_calc_ratio_first_second_add_mean() {
        // The scaled signals themselves, without pseudocounts.
        let cov1: f32 = 5.0;
        let cov2: f32 = 10.0;
        let sf1: f32 = 1.0;
        let sf2: f32 = 2.0;
        let pc1: f32 = 1.0;
        let pc2: f32 = 2.0;
        assert_eq!(calc_ratio(cov1, cov2, &sf1, &sf2, &pc1, &pc2, "first"), 5.0);
        assert_eq!(calc_ratio(cov1, cov2, &sf1, &sf2, &pc1, &pc2, "second"), 20.0);
        assert_eq!(calc_ratio(cov1, cov2, &sf1, &sf2, &pc1, &pc2, "add"), 25.0);
        assert_eq!(calc_ratio(cov1, cov2, &sf1, &sf2, &pc1, &pc2, "mean"), 12.5);
        // and none of them is the log2 ratio
        let l2 = calc_ratio(cov1, cov2, &sf1, &sf2, &pc1, &pc2, "log2");
        for op in ["first", "second", "add", "mean"] {
            assert_ne!(calc_ratio(cov1, cov2, &sf1, &sf2, &pc1, &pc2, op), l2);
        }
    }

    #[test]
    fn test_deseq_scalefactors() {
        let counts =
            Array2::from_shape_vec((3, 3), vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
                .unwrap();
        let sf = deseq_scalefactors(&counts);
        let expected_sf = Array1::from_vec(vec![1.233106, 0.9864848, 0.82207066]);
        assert_eq!(sf, expected_sf);
    }
}
