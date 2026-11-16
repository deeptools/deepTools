//! Parallel reader for computeMatrix's gzipped text matrix (used by
//! `heatmapper.read_matrix_file`, i.e. plotHeatmap / plotProfile /
//! computeMatrixOperations).
//!
//! Semantics follow the per-line Python reader: each line is stripped, lines
//! starting with '@' carry the JSON parameters, other lines are split on tabs
//! into 6 region fields and the values. Values are parsed with
//! `str::parse::<f64>`, which is correctly rounded like Python's `float()`
//! (and accepts nan/inf). Anything the fast path does not handle (a parameter
//! line after data rows, ragged rows, unparsable tokens) is reported as a
//! ValueError so the caller can fall back to the Python reader.

use flate2::read::MultiGzDecoder;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyByteArray;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::mpsc::sync_channel;

/// Decompressed bytes per parse chunk.
const CHUNK_BYTES: usize = 32 << 20;

type Meta = (String, String, String, String, String, String);

fn parse_line(line: &[u8]) -> Result<(Meta, Vec<f64>), String> {
    let line = std::str::from_utf8(line)
        .map_err(|_| "non UTF-8 line".to_string())?
        .trim_matches(|c: char| c.is_ascii_whitespace() || c == '\x0b');
    let mut it = line.split('\t');
    let mut f = || {
        it.next()
            .map(|s| s.to_string())
            .ok_or_else(|| "line with fewer than 6 fields".to_string())
    };
    let meta = (f()?, f()?, f()?, f()?, f()?, f()?);
    let vals = it
        .map(|t| {
            t.trim_ascii()
                .parse::<f64>()
                .map_err(|_| format!("cannot parse value '{}'", t))
        })
        .collect::<Result<Vec<f64>, String>>()?;
    Ok((meta, vals))
}

/// Returns (parameter_lines, region_fields, values_as_f64_bytes, nrows, ncols).
#[pyfunction]
pub fn r_read_matrix<'py>(
    py: Python<'py>,
    path: &str,
    nproc: usize,
) -> PyResult<(
    Vec<String>,
    Vec<Meta>,
    Bound<'py, PyByteArray>,
    usize,
    usize,
)> {
    let path = path.to_string();
    let res: Result<(Vec<String>, Vec<Meta>, Vec<f64>, usize), String> = py.detach(move || {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(nproc.max(1))
            .build()
            .map_err(|e| e.to_string())?;
        let file = File::open(&path).map_err(|e| format!("{}: {}", path, e))?;
        // Decompress on a separate thread, handing over chunks of whole lines.
        let (tx, rx) = sync_channel::<Result<Vec<Vec<u8>>, String>>(2);
        let reader = std::thread::spawn(move || {
            let mut rd = BufReader::with_capacity(1 << 20, MultiGzDecoder::new(file));
            loop {
                let mut lines: Vec<Vec<u8>> = Vec::new();
                let mut bytes = 0usize;
                while bytes < CHUNK_BYTES {
                    let mut l = Vec::new();
                    match rd.read_until(b'\n', &mut l) {
                        Ok(0) => break,
                        Ok(n) => {
                            bytes += n;
                            lines.push(l);
                        }
                        Err(e) => {
                            let _ = tx.send(Err(e.to_string()));
                            return;
                        }
                    }
                }
                let done = lines.is_empty();
                if !done && tx.send(Ok(lines)).is_err() {
                    return;
                }
                if done || bytes < CHUNK_BYTES {
                    return;
                }
            }
        });
        let mut params: Vec<String> = Vec::new();
        let mut metas: Vec<Meta> = Vec::new();
        let mut values: Vec<f64> = Vec::new();
        let mut ncols: Option<usize> = None;
        for chunk in rx {
            let mut lines = chunk?;
            // Parameter lines ('@') are only allowed before the first data row.
            let mut skip = 0;
            while skip < lines.len() && lines[skip].trim_ascii_start().first() == Some(&b'@') {
                if !metas.is_empty() {
                    return Err("parameter line after data rows".to_string());
                }
                let l = String::from_utf8(std::mem::take(&mut lines[skip]))
                    .map_err(|_| "non UTF-8 parameter line".to_string())?;
                params.push(l);
                skip += 1;
            }
            let parsed: Vec<(Meta, Vec<f64>)> = pool.install(|| {
                lines[skip..]
                    .par_iter()
                    .map(|l| {
                        if l.trim_ascii_start().first() == Some(&b'@') {
                            Err("parameter line after data rows".to_string())
                        } else {
                            parse_line(l)
                        }
                    })
                    .collect::<Result<Vec<_>, String>>()
            })?;
            for (meta, row) in parsed {
                match ncols {
                    None => ncols = Some(row.len()),
                    Some(n) if n != row.len() => return Err("rows of unequal length".to_string()),
                    _ => {}
                }
                values.extend_from_slice(&row);
                metas.push(meta);
            }
        }
        reader
            .join()
            .map_err(|_| "decompression thread panicked".to_string())?;
        Ok((params, metas, values, ncols.unwrap_or(0)))
    });
    let (params, metas, values, ncols) = res.map_err(PyValueError::new_err)?;
    let nrows = metas.len();
    let nbytes = values.len() * std::mem::size_of::<f64>();
    let buf = PyByteArray::new_with(py, nbytes, |b: &mut [u8]| {
        for (dst, v) in b.chunks_exact_mut(8).zip(values.iter()) {
            dst.copy_from_slice(&v.to_ne_bytes());
        }
        Ok(())
    })?;
    Ok((params, metas, buf, nrows, ncols))
}
