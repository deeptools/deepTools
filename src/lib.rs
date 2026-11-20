//! deepTools core. The `python` feature (default) builds the Python extension.
#![warn(missing_docs)]

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg_attr(not(feature = "python"), allow(dead_code))]
mod calc;
/// Region parsing and per-bin BAM pileup.
#[cfg(feature = "io")]
#[allow(missing_docs)]
pub mod covcalc;
/// BAM, BED, GTF and bigWig readers and writers.
#[cfg(feature = "io")]
#[allow(missing_docs)]
pub mod filehandler;
/// Per-alignment filters (MAPQ, flags, fragment length, blacklist, extension).
#[cfg(feature = "io")]
#[allow(missing_docs)]
pub mod filtering;
#[cfg(feature = "python")]
mod matrix;
pub mod normalization;

pub use calc::deseq_scalefactors;

/// Fragment shifting of BAM records (`apply_shift`). The python command lives here too.
#[cfg(feature = "io")]
#[allow(missing_docs)]
pub mod alignmentsieve;
#[cfg(feature = "python")]
mod bamcompare;
#[cfg(feature = "python")]
mod bamcoverage;
#[cfg(feature = "python")]
mod computematrix;
#[cfg(feature = "python")]
mod multibamsummary;

#[cfg(test)]
mod tests;

#[cfg(feature = "python")]
#[pymodule]
fn hp(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(bamcoverage::r_bamcoverage, m)?)?;
    m.add_function(wrap_pyfunction!(bamcompare::r_bamcompare, m)?)?;
    m.add_function(wrap_pyfunction!(computematrix::r_computematrix, m)?)?;
    m.add_function(wrap_pyfunction!(alignmentsieve::r_alignmentsieve, m)?)?;
    m.add_function(wrap_pyfunction!(multibamsummary::r_mbams, m)?)?;
    m.add_function(wrap_pyfunction!(matrix::reader::r_read_matrix, m)?)?;
    Ok(())
}
