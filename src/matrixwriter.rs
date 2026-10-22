//! Text formatting of computeMatrix's matrix rows.
//!
//! Values are formatted with `fixedfmt::push_f32_fixed6`, which reproduces
//! `format!("{:.6}", v)` exactly, into a reusable buffer instead of allocating a
//! `String` per value.

use crate::covcalc::Region;
use crate::fixedfmt::push_f32_fixed6;
use std::io::Write;

/// Append the 6 BED-like fields and the trailing tab.
#[inline]
pub fn push_region_fields(buf: &mut Vec<u8>, region: &Region) {
    write!(
        buf,
        "{}\t{}\t{}\t{}\t{}\t{}\t",
        region.chrom,            // Chromosome
        region.start.rewrites(), // Revalue for start (either u32, or Vec<u32>)
        region.end.rewrites(),   // Revalue for end (either u32, or Vec<u32>)
        region.name,  // String field for name. Duplicates taken care of in read_bedfiles.
        region.score, // Score field persisted from bedfile
        region.strand, // Strand field persisted from bedfile
    )
    .expect("writing to a Vec cannot fail");
}

/// Append the tab-separated values of one row, without a trailing newline.
#[inline]
pub fn push_values(buf: &mut Vec<u8>, row: &[f32], scale: f32) {
    for (i, x) in row.iter().enumerate() {
        if i > 0 {
            buf.push(b'\t');
        }
        if x.is_nan() {
            // as deepTools 3 encoded nan in matrix.
            buf.extend_from_slice(b"nan");
        } else {
            push_f32_fixed6(buf, scale * x);
        }
    }
}

/// Append one full matrix line.
#[inline]
pub fn push_row(buf: &mut Vec<u8>, region: &Region, row: &[f32], scale: f32) {
    push_region_fields(buf, region);
    push_values(buf, row, scale);
    buf.push(b'\n');
}
