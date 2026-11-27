//! deepTools core. The `python` feature (default) builds the Python extension.
#![warn(missing_docs)]

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg_attr(not(feature = "python"), allow(dead_code))]
mod calc;
#[cfg(feature = "io")]
#[cfg_attr(not(feature = "python"), allow(dead_code))]
mod covcalc;
#[cfg(feature = "io")]
#[cfg_attr(not(feature = "python"), allow(dead_code))]
mod filehandler;
#[cfg(feature = "io")]
#[cfg_attr(not(feature = "python"), allow(dead_code))]
mod filtering;
#[cfg(feature = "python")]
mod matrix;
pub mod normalization;

pub use calc::deseq_scalefactors;

/// BAM, BED and bigWig I/O for Rust programs that embed deepTools without the Python extension
/// (build with `--no-default-features --features io`).
///
/// This is a deliberately small surface: only what such programs need to read regions, filter and shift
/// alignments and write coverage files. Everything else in the crate stays private. Items here are
/// not covered by any stability promise and may change in any release.
#[cfg(feature = "io")]
pub mod io {
    pub use crate::alignmentsieve::{ShiftRefusal, apply_shift};
    pub use crate::filehandler::{read_bedfile, write_covfile};
    pub use crate::filtering::Alignmentfilters;
}

#[cfg(feature = "io")]
#[cfg_attr(not(feature = "python"), allow(dead_code))]
mod alignmentsieve;
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
