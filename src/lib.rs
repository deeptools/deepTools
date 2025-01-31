use pyo3::prelude::*;
mod alignmentsieve;
mod bamcoverage;
mod bamcompare;
mod calc;
mod computematrix;
mod covcalc;
mod filehandler;
mod multibamsummary;
mod normalization;
mod test;

#[pymodule]
fn hp(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(bamcoverage::r_bamcoverage, m)?)?;
    m.add_function(wrap_pyfunction!(bamcompare::r_bamcompare, m)?)?;
    m.add_function(wrap_pyfunction!(computematrix::r_computematrix, m)?)?;
    m.add_function(wrap_pyfunction!(alignmentsieve::r_alignmentsieve, m)?)?;
    m.add_function(wrap_pyfunction!(multibamsummary::r_mbams, m)?)?;
    Ok(())
}
