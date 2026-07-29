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
        assert_eq!(r3, -0.27);
        assert_eq!(r4, -16.0);
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
