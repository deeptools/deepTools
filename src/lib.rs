#[cfg(feature = "python")]
use pyo3::prelude::*;
mod calc;
mod covcalc;
mod filehandler;
mod filtering;
mod normalization;

#[cfg(feature = "python")]
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
    Ok(())
}