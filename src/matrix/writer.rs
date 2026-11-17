//! computeMatrix output: the gzipped matrix, `--outFileNameMatrix` and `--outFileSortedRegions`.
//!
//! The gzipped matrix is written in parallel. Values are formatted with
//! `fixedfmt::push_f32_fixed6`, which reproduces `format!("{:.6}", v)` exactly,
//! and the text is compressed pigz-style:
//!
//! * rows are grouped into ~1 MiB text chunks; a batch of chunks is formatted in
//!   parallel and then deflated in parallel (rayon, current pool);
//! * each chunk is raw-deflated independently and ended with a sync flush
//!   (the last one with a final block), primed with the previous chunk's last
//!   32 KiB as a preset dictionary, so the concatenation is one valid deflate
//!   stream inside a single gzip member (plain `gzip`/`zcat`/Python `gzip.open`
//!   read it, and `gzip -l` reports the true size);
//! * the CRC32 of the whole text is assembled with `crc32fast::Hasher::combine`.

use super::fixedfmt::push_f32_fixed6;
use crate::covcalc::{Region, Revalue, Scalingregions};
use flate2::{Compress, Compression, FlushCompress, Status};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

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

pub fn write_matrix(
    header: String,
    mat: &[Vec<f32>],
    ofile: &str,
    regions: &[Region],
    scale_regions: &Scalingregions,
) {
    // Write out the matrix to a compressed file. Rows are formatted and
    // compressed in parallel on the current rayon pool (write_matrix_gz).
    let omat = File::create(ofile)
        .unwrap_or_else(|e| panic!("Failed to create output matrix file '{}': {}", ofile, e));
    write_matrix_gz(
        &mut BufWriter::with_capacity(1 << 20, omat),
        header.as_bytes(),
        mat,
        regions,
        scale_regions.scale,
        Compression::default(),
    )
    .unwrap_or_else(|e| panic!("Failed to write matrix file '{}': {}", ofile, e));
}

pub fn write_matrix_values(
    file_name: &str,
    mat: &[Vec<f32>],
    scale_regions: &Scalingregions,
    regionsizes: &HashMap<String, u32>,
) {
    use std::io::Write;
    let mut fh = File::create(file_name)
        .unwrap_or_else(|e| panic!("Failed to create matrix values file '{}': {}", file_name, e));

    // Header line 1: group labels with region counts
    let info: Vec<String> = scale_regions
        .regionlabels
        .iter()
        .map(|label| format!("{}:{}", label, regionsizes.get(label).unwrap_or(&0)))
        .collect();
    fh.write_all(format!("#{}\n", info.join("\t")).as_bytes())
        .unwrap_or_else(|e| {
            panic!(
                "Failed to write header to matrix values file '{}': {}",
                file_name, e
            )
        });

    // Header line 2: region dimension parameters
    let header2 = format!(
        "#downstream:{}\tupstream:{}\tbody:{}\tbin size:{}\tunscaled 5 prime:{}\tunscaled 3 prime:{}\n",
        scale_regions.downstream,
        scale_regions.upstream,
        scale_regions.regionbodylength,
        scale_regions.binsize,
        scale_regions.unscaled5prime,
        scale_regions.unscaled3prime,
    );
    fh.write_all(header2.as_bytes()).unwrap_or_else(|e| {
        panic!(
            "Failed to write header line 2 to matrix values file '{}': {}",
            file_name, e
        )
    });

    // Header line 3: sample labels repeated per column
    let cols_per_sample = scale_regions.cols_expected / scale_regions.bwfiles;
    let sample_info: Vec<String> = scale_regions
        .bwlabels
        .iter()
        .flat_map(|label| (0..cols_per_sample).map(move |_| label.clone()))
        .collect();
    fh.write_all(format!("{}\t{}\n", info.join("\t"), sample_info.join("\t")).as_bytes())
        .unwrap_or_else(|e| {
            panic!(
                "Failed to write header line 3 to matrix values file '{}': {}",
                file_name, e
            )
        });
    fh.flush()
        .unwrap_or_else(|e| panic!("Failed to flush matrix values file '{}': {}", file_name, e));

    // Reopen in append mode and write matrix data
    let fh = std::fs::OpenOptions::new()
        .create(false)
        .append(true)
        .open(file_name)
        .unwrap_or_else(|e| {
            panic!(
                "Failed to reopen matrix values file '{}' for appending: {}",
                file_name, e
            )
        });
    let mut fh = BufWriter::with_capacity(1 << 20, fh);
    let mut line: Vec<u8> = Vec::new();
    for row in mat.iter() {
        line.clear();
        push_values(&mut line, row, scale_regions.scale);
        line.push(b'\n');
        fh.write_all(&line)
            .unwrap_or_else(|e| panic!("Failed to write matrix row to '{}': {}", file_name, e));
    }
    fh.flush()
        .unwrap_or_else(|e| panic!("Failed to flush matrix values file '{}': {}", file_name, e));
}

pub fn write_sorted_regions_bed(
    file_name: &str,
    regions: &[Region],
    scale_regions: &Scalingregions,
    regionsizes: &HashMap<String, u32>,
) {
    use std::io::Write;
    let mut fh = File::create(file_name).unwrap_or_else(|e| {
        panic!(
            "Failed to create sorted regions BED file '{}': {}",
            file_name, e
        )
    });
    let header = "#chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRGB\tblockCount\tblockSizes\tblockStarts\tdeepTools_group\n";
    fh.write_all(header.as_bytes()).unwrap_or_else(|e| {
        panic!(
            "Failed to write header to sorted regions BED file '{}': {}",
            file_name, e
        )
    });
    // Build group boundaries from regionsizes (cumulative region counts per group label)
    let mut group_boundaries: Vec<u32> = vec![0];
    let mut cumsum: u32 = 0;
    for label in &scale_regions.regionlabels {
        cumsum += *regionsizes.get(label).unwrap_or(&0);
        group_boundaries.push(cumsum);
    }
    for (idx, region) in regions.iter().enumerate() {
        let label_idx = group_boundaries
            .iter()
            .take_while(|&&b| b <= idx as u32)
            .count()
            .saturating_sub(1)
            .min(scale_regions.regionlabels.len() - 1);
        let start_first = match &region.start {
            Revalue::U(v) => *v,
            Revalue::V(vs) => *vs.first().unwrap_or_else(|| {
                panic!("Region '{}' has an empty exon-start vector", region.name)
            }),
        };
        let end_last = match &region.end {
            Revalue::U(v) => *v,
            Revalue::V(vs) => *vs
                .last()
                .unwrap_or_else(|| panic!("Region '{}' has an empty exon-end vector", region.name)),
        };
        let (block_count, block_sizes, block_starts) = match (&region.start, &region.end) {
            (Revalue::U(_), _) => (
                "1".to_string(),
                (end_last - start_first).to_string(),
                "0".to_string(),
            ),
            (Revalue::V(starts), Revalue::V(ends)) => {
                let bc = starts.len().to_string();
                let sz: Vec<String> = starts
                    .iter()
                    .zip(ends.iter())
                    .map(|(s, e)| (e - s).to_string())
                    .collect();
                let st: Vec<String> = starts
                    .iter()
                    .map(|s| (s - start_first).to_string())
                    .collect();
                (bc, sz.join(","), st.join(","))
            }
            _ => (
                "1".to_string(),
                (end_last - start_first).to_string(),
                "0".to_string(),
            ),
        };
        let group_label = &scale_regions.regionlabels[label_idx];
        let score_str = region
            .score
            .parse::<f64>()
            .map(|v| format!("{:.1}", v))
            .unwrap_or_else(|_| region.score.clone());

        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t{}\n",
            region.chrom,
            start_first,
            end_last,
            region.name,
            score_str,
            region.strand,
            start_first,
            end_last,
            block_count,
            block_sizes,
            block_starts,
            group_label,
        );
        fh.write_all(line.as_bytes()).unwrap_or_else(|e| {
            panic!(
                "Failed to write line to sorted regions BED file '{}': {}",
                file_name, e
            )
        });
    }
}
