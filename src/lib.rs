//! deepTools core. The `python` feature (default) builds the Python extension.
#![warn(missing_docs)]

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg_attr(not(feature = "python"), allow(dead_code))]
mod calc;
#[cfg(feature = "python")]
mod covcalc;
#[cfg(feature = "python")]
mod filehandler;
#[cfg(feature = "python")]
mod filtering;
#[cfg(feature = "python")]
mod fixedfmt;
#[cfg(feature = "python")]
mod matrixwriter;
pub mod normalization;

pub use calc::deseq_scalefactors;

#[cfg(feature = "python")]
mod alignmentsieve;
#[cfg(feature = "python")]
mod bamcompare;
#[cfg(feature = "python")]
mod bamcoverage;
#[cfg(feature = "python")]
mod computematrix;
#[cfg(feature = "python")]
mod matrixreader;
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
    m.add_function(wrap_pyfunction!(matrixreader::r_read_matrix, m)?)?;
    Ok(())
}
