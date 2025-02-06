use ndarray::{Array1, Array2, Axis};

pub fn median(mut nvec: Vec<u32>) -> f32 {
    if nvec.is_empty() {
        0.0
    } else if nvec.len() == 1 {
        nvec[0] as f32
    } else {
        nvec.sort_unstable();
        let len = nvec.len();
        if len % 2 == 1 {
            nvec[len / 2] as f32
        } else {
            (nvec[len / 2] + nvec[len / 2 - 1]) as f32 / 2.0
        }
    }
}

pub fn mean_float(fvec: Vec<&f32>) -> f32 {
    let valid_floats: Vec<f32> = fvec
        .into_iter()
        .cloned()
        .filter(|v| v.is_finite())
        .collect();
    if valid_floats.is_empty() {
        0.0
    } else {
        valid_floats.iter().sum::<f32>() / valid_floats.len() as f32
    }
}

pub fn median_float(fvec: Vec<&f32>) -> f32 {
    let mut valid_floats: Vec<f32> = fvec
        .into_iter()
        .cloned()
        .filter(|v| v.is_finite())
        .collect();
    if valid_floats.is_empty() {
        0.0
    } else {
        valid_floats.sort_by(
            |a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
        );
        let len = valid_floats.len();
        if len % 2 == 1 {
            valid_floats[len / 2]
        } else {
            (valid_floats[len / 2] + valid_floats[len / 2 - 1]) / 2.0
        }
    }
}

pub fn min_float(fvec: Vec<&f32>) -> f32 {
    let valid_floats: Vec<f32> = fvec
        .into_iter()
        .cloned()
        .filter(|v| v.is_finite())
        .collect();
    if valid_floats.is_empty() {
        0.0
    } else {
        valid_floats
            .iter()
            .cloned()
            .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Greater))
            .unwrap_or(0.0)
    }
}

pub fn max_float(fvec: Vec<&f32>) -> f32 {
    let valid_floats: Vec<f32> = fvec
        .into_iter()
        .cloned()
        .filter(|v| v.is_finite())
        .collect();
    if valid_floats.is_empty() {
        0.0
    } else {
        valid_floats
            .iter()
            .cloned()
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Greater))
            .unwrap_or(0.0)
    }
}

pub fn sum_float(fvec: Vec<&f32>) -> f32 {
    let valid_floats: Vec<f32> = fvec
        .into_iter()
        .cloned()
        .filter(|v| v.is_finite())
        .collect();
    if valid_floats.is_empty() {
        0.0
    } else {
        valid_floats.iter().copied().sum::<f32>()
    }
}

pub fn std_float(fvec: Vec<&f32>) -> f32 {
    let valid_floats: Vec<f32> = fvec
        .into_iter()
        .cloned()
        .filter(|v| v.is_finite())
        .collect();
    if valid_floats.is_empty() {
        0.0
    } else {
        let n = valid_floats.len() as f32;
        let mean = valid_floats.iter().copied().sum::<f32>() / n;
        let stdsum = valid_floats
            .iter()
            .copied()
            .map(|val| (val - mean).powi(2) as f64)  // Squared difference from mean
            .sum::<f64>();
        let stdsumnorm = stdsum / (n - 1.0) as f64;
        let stdsumnormsqrt = stdsumnorm.sqrt();
        stdsumnormsqrt as f32
    }
}

pub fn calc_ratio(
    cov1: f32,
    cov2: f32,
    sf1: &f32,
    sf2: &f32,
    pseudocount: &f32,
    operation: &str
) -> f32 {
    // Pseudocounts are only used in log2 and ratio operations
    // First scale factor is applied, then pseudocount, if applicable.
    match operation {
        "log2" => {
            let num: f32 = (cov1 * *sf1) + *pseudocount;
            let den: f32 = (cov2 * *sf2) + *pseudocount;
            return (num / den).log2();
        }
        "ratio" => {
            let num: f32 = (cov1 * *sf1) + *pseudocount;
            let den: f32 = (cov2 * *sf2) + *pseudocount;
            return num / den;
        }
        "reciprocal_ratio" => {
            let num: f32 = (cov1 * *sf1) + *pseudocount;
            let den: f32 = (cov2 * *sf2) + *pseudocount;
            let ratio: f32 = num / den;
            if ratio >= 1.0 {
                return den / num;
            } else {
                return -num / den;
            }
        }
        _ => {
            // No operation is never allowed (on the py arg level, so just default to log2)
            let num: f32 = (cov1 * *sf1) + *pseudocount;
            let den: f32 = (cov2 * *sf2) + *pseudocount;
            return (num / den).log2();
        }
    }
}

pub fn deseq_scalefactors(array2: &Array2<f32>) -> Array1<f32> {
    let loggeomeans = array2.mapv(|v| v.ln()).sum_axis(Axis(1)) / array2.shape()[1] as f32;
    let masked_array = array2.mapv(|x| if x <= 0.0 { f32::NAN } else { x });
    let masked_loggeomeans = loggeomeans.mapv(|x| if x.is_infinite() { f32::NAN } else { x });
    let adjusted_loga = masked_array.mapv(|x| x.ln()).t().to_owned() - &masked_loggeomeans;
    let medians: Array1<f32> = adjusted_loga.t().axis_iter(Axis(1))
        .map(|x| {
            let vec: Vec<&f32> = x.iter()
                .filter(|&&x| !x.is_nan())
                .collect();
            median_float(vec)
        })
        .collect();
    medians.mapv(|x| 1.0 / x.exp())
}