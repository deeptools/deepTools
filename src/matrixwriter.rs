//! Parallel writer for computeMatrix's gzipped text matrix.
//!
//! Values are formatted with `fixedfmt::push_f32_fixed6`, which reproduces
//! `format!("{:.6}", v)` exactly, and the text is compressed pigz-style:
//!
//! * rows are grouped into ~1 MiB text chunks; a batch of chunks is formatted in
//!   parallel and then deflated in parallel (rayon, current pool);
//! * each chunk is raw-deflated independently and ended with a sync flush
//!   (the last one with a final block), primed with the previous chunk's last
//!   32 KiB as a preset dictionary, so the concatenation is one valid deflate
//!   stream inside a single gzip member (plain `gzip`/`zcat`/Python `gzip.open`
//!   read it, and `gzip -l` reports the true size);
//! * the CRC32 of the whole text is assembled with `crc32fast::Hasher::combine`.

use crate::covcalc::Region;
use crate::fixedfmt::push_f32_fixed6;
use flate2::{Compress, Compression, FlushCompress, Status};
use rayon::prelude::*;
use std::io::Write;

/// Uncompressed bytes per independently compressed chunk.
const CHUNK_BYTES: usize = 1 << 20;
/// deflate window size: the preset dictionary length used to prime each chunk.
const DICT_BYTES: usize = 32 * 1024;

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

/// Raw-deflate `input` as one piece of a larger stream. Non-final pieces end with a
/// sync flush (byte aligned, no BFINAL), so pieces can be concatenated.
fn deflate_piece(input: &[u8], dict: &[u8], level: Compression, last: bool) -> Vec<u8> {
    let mut c = Compress::new(level, false);
    if !dict.is_empty() {
        c.set_dictionary(dict)
            .expect("failed to set deflate dictionary");
    }
    let flush = if last {
        FlushCompress::Finish
    } else {
        FlushCompress::Sync
    };
    let mut out: Vec<u8> = Vec::with_capacity(input.len() / 4 + 1024);
    loop {
        if out.capacity() - out.len() < 1024 {
            out.reserve(out.capacity() / 2 + 4096);
        }
        let consumed = c.total_in() as usize;
        let status = c
            .compress_vec(&input[consumed..], &mut out, flush)
            .expect("deflate failed");
        let all_in = c.total_in() as usize == input.len();
        match status {
            Status::StreamEnd => break,
            // zlib semantics: a flush is complete once the call returns with
            // output space left over.
            _ if !last && all_in && out.len() < out.capacity() => break,
            _ => {}
        }
    }
    out
}

/// Gzip member header, identical to what `flate2::write::GzEncoder` writes for
/// `Compression::default()` (mtime 0, XFL 0, OS 255).
const GZ_HEADER: [u8; 10] = [0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 255];

/// Write `header` followed by one text line per (region, row) as a single gzip
/// member, formatting and compressing in parallel on the current rayon pool.
pub fn write_matrix_gz<W: Write>(
    out: &mut W,
    header: &[u8],
    mat: &[Vec<f32>],
    regions: &[Region],
    scale: f32,
    level: Compression,
) -> std::io::Result<()> {
    assert_eq!(regions.len(), mat.len());
    let nrows = mat.len();
    let ncols = mat.first().map_or(0, |r| r.len());
    let est_row_bytes = 64 + 10 * ncols;
    let rows_per_chunk = (CHUNK_BYTES / est_row_bytes).max(1);
    let row_chunks: Vec<(usize, usize)> = (0..nrows)
        .step_by(rows_per_chunk)
        .map(|s| (s, (s + rows_per_chunk).min(nrows)))
        .collect();
    // Bounds the memory held at once to a few chunks per thread.
    let batch = (4 * rayon::current_num_threads()).max(1);

    out.write_all(&GZ_HEADER)?;
    let mut crc = crc32fast::Hasher::new();
    let mut total_len: u64 = 0;

    // The JSON header is the first piece of the stream.
    out.write_all(&deflate_piece(header, &[], level, row_chunks.is_empty()))?;
    crc.update(header);
    total_len += header.len() as u64;
    let mut prev_tail: Vec<u8> = header[header.len().saturating_sub(DICT_BYTES)..].to_vec();

    let nchunks = row_chunks.len();
    for (bi, batch_chunks) in row_chunks.chunks(batch).enumerate() {
        let texts: Vec<Vec<u8>> = batch_chunks
            .par_iter()
            .map(|&(s, e)| {
                let mut buf = Vec::with_capacity((e - s) * est_row_bytes + 256);
                for i in s..e {
                    push_row(&mut buf, &regions[i], &mat[i], scale);
                }
                buf
            })
            .collect();
        let first_idx = bi * batch;
        let compressed: Vec<(Vec<u8>, crc32fast::Hasher)> = texts
            .par_iter()
            .enumerate()
            .map(|(j, text)| {
                let dict: &[u8] = if j == 0 {
                    &prev_tail
                } else {
                    let p = &texts[j - 1];
                    &p[p.len().saturating_sub(DICT_BYTES)..]
                };
                let last = first_idx + j + 1 == nchunks;
                let mut h = crc32fast::Hasher::new();
                h.update(text);
                (deflate_piece(text, dict, level, last), h)
            })
            .collect();
        for ((bytes, h), text) in compressed.into_iter().zip(texts.iter()) {
            out.write_all(&bytes)?;
            crc.combine(&h);
            total_len += text.len() as u64;
        }
        if let Some(t) = texts.last() {
            prev_tail = t[t.len().saturating_sub(DICT_BYTES)..].to_vec();
        }
    }
    // Gzip trailer: CRC32 and the uncompressed size modulo 2^32.
    out.write_all(&crc.finalize().to_le_bytes())?;
    out.write_all(&(total_len as u32).to_le_bytes())?;
    out.flush()
}
