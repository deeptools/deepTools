use ndarray::{Array1, Array2, Axis};

pub fn median(mut nvec: Vec<u32>) -> f32 {
    if nvec.is_empty() {
        0.0
    } else if nvec.len() == 1 {
        nvec[0] as f32
    } else {
        let len = nvec.len();
        if len % 2 == 1 {
            let mid = len / 2;
            nvec.select_nth_unstable(mid);
            nvec[mid] as f32
        } else {
            let m = len / 2;
            nvec.select_nth_unstable(m - 1);
            nvec[m..].select_nth_unstable(0);
            (nvec[m - 1] + nvec[m]) as f32 / 2.0
        }
    }
}

pub fn mean_float(fvec: &[&f32]) -> f32 {
    let mut sum = 0.0f64;
    let mut count = 0usize;
    for &&v in fvec {
        if v.is_finite() {
            sum += v as f64;
            count += 1;
        }
    }
    if count == 0 { 0.0 } else { (sum / count as f64) as f32 }
}

pub fn median_float(fvec: &[&f32]) -> f32 {
    let mut valid_floats: Vec<f32> = fvec
        .iter()
        .copied()
        .copied()
        .filter(|v| v.is_finite())
        .collect();
    if valid_floats.is_empty() {
        return 0.0;
    }
    valid_floats.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let len = valid_floats.len();
    if len % 2 == 1 {
        valid_floats[len / 2]
    } else {
        (valid_floats[len / 2] + valid_floats[len / 2 - 1]) / 2.0
    }
}

pub fn min_float(fvec: &[&f32]) -> f32 {
    let mut min_val = f32::NAN;
    for &&v in fvec {
        if v.is_finite() {
            if min_val.is_nan() || v < min_val {
                min_val = v;
            }
        }
    }
    if min_val.is_nan() { 0.0 } else { min_val }
}

pub fn max_float(fvec: &[&f32]) -> f32 {
    let mut max_val = f32::NAN;
    for &&v in fvec {
        if v.is_finite() {
            if max_val.is_nan() || v > max_val {
                max_val = v;
            }
        }
    }
    if max_val.is_nan() { 0.0 } else { max_val }
}

pub fn sum_float(fvec: &[&f32]) -> f32 {
    let mut sum = 0.0f32;
    let mut count = false;
    for &&v in fvec {
        if v.is_finite() {
            sum += v;
            count = true;
        }
    }
    if count { sum } else { 0.0 }
}

pub fn std_float(fvec: &[&f32]) -> f32 {
    let valid_floats: Vec<f64> = fvec
        .iter()
        .copied()
        .copied()
        .filter(|v| v.is_finite())
        .map(|v| v as f64)
        .collect();
    if valid_floats.is_empty() {
        return 0.0f32;
    }
    let n = valid_floats.len() as f64;
    let mean = valid_floats.iter().copied().sum::<f64>() / n;
    let stdsum = valid_floats
        .iter()
        .copied()
        .map(|val| (val - mean).powi(2))
        .sum::<f64>();
    let stdsumnorm = stdsum / n; // population std
    stdsumnorm.sqrt() as f32
}

pub fn calc_ratio(
    cov1: f32,
    cov2: f32,
    sf1: &f32,
    sf2: &f32,
    pseudocount1: &f32,
    pseudocount2: &f32,
    operation: &str,
) -> f32 {
    // Pseudocounts are only used in log2 and ratio operations
    // First scale factor is applied, then pseudocount, if applicable.
    match operation {
        "log2" => {
            let num: f32 = (cov1 * *sf1) + *pseudocount1;
            let den: f32 = (cov2 * *sf2) + *pseudocount2;
            let fcov: f32 = (num / den).log2();
            return (fcov * 100.0).round() / 100.0;
        }
        "ratio" => {
            let num: f32 = (cov1 * *sf1) + *pseudocount1;
            let den: f32 = (cov2 * *sf2) + *pseudocount2;
            let fcov: f32 = num / den;
            return (fcov * 100.0).round() / 100.0;
        }
        "reciprocal_ratio" => {
            // a/b if a/b >= 1, else -b/a (negative fold change), as in the
            // --operation help and the 3.5.x implementation.
            let num: f32 = (cov1 * *sf1) + *pseudocount1;
            let den: f32 = (cov2 * *sf2) + *pseudocount2;
            let ratio: f32 = num / den;
            if ratio >= 1.0 {
                return (ratio * 100.0).round() / 100.0;
            } else {
                let fcov: f32 = -den / num;
                return (fcov * 100.0).round() / 100.0;
            }
        }
        "subtract" => {
            let num: f32 = (cov1 * *sf1) + *pseudocount1;
            let den: f32 = (cov2 * *sf2) + *pseudocount2;
            let fcov: f32 = num - den;
            return (fcov * 100.0).round() / 100.0;
        }
        // The remaining operations output the scaled signal(s) without any
        // pseudocount (the --pseudocount help: only used with log2 / ratio).
        "first" => {
            let fcov: f32 = cov1 * *sf1;
            return (fcov * 100.0).round() / 100.0;
        }
        "second" => {
            let fcov: f32 = cov2 * *sf2;
            return (fcov * 100.0).round() / 100.0;
        }
        "add" => {
            let fcov: f32 = (cov1 * *sf1) + (cov2 * *sf2);
            return (fcov * 100.0).round() / 100.0;
        }
        "mean" => {
            let fcov: f32 = ((cov1 * *sf1) + (cov2 * *sf2)) / 2.0;
            return (fcov * 100.0).round() / 100.0;
        }
        _ => {
            // The CLI restricts --operation to the eight choices above.
            panic!("Unknown bamCompare operation '{}'", operation);
        }
    }
}

pub fn deseq_scalefactors(array2: &Array2<f32>) -> Array1<f32> {
    let loggeomeans = array2.mapv(|v| v.ln()).sum_axis(Axis(1)) / array2.shape()[1] as f32;
    let masked_array = array2.mapv(|x| if x <= 0.0 { f32::NAN } else { x });
    let masked_loggeomeans = loggeomeans.mapv(|x| if x.is_infinite() { f32::NAN } else { x });
    let adjusted_loga = masked_array.mapv(|x| x.ln()).t().to_owned() - &masked_loggeomeans;
    let medians: Array1<f32> = adjusted_loga
        .t()
        .axis_iter(Axis(1))
        .map(|x| {
            let vec: Vec<&f32> = x.iter().filter(|&&x| !x.is_nan()).collect();
            median_float(&vec)
        })
        .collect();
    medians.mapv(|x| 1.0 / x.exp())
}
