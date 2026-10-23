use crate::covcalc::{Region, Revalue};
use crate::matrixwriter::write_matrix_gz;
use crate::tests::test_fixedfmt::Rng;
use flate2::Compression;
use flate2::read::GzDecoder;
use std::io::Read;

fn test_regions_matrix(nrows: usize, ncols: usize, seed: u64) -> (Vec<Region>, Vec<Vec<f32>>) {
    let mut rng = Rng(seed);
    let mut regions = Vec::with_capacity(nrows);
    let mut mat = Vec::with_capacity(nrows);
    for i in 0..nrows {
        let start = (rng.next() % 1_000_000) as u32;
        let (s, e) = if i % 3 == 0 {
            (
                Revalue::V(vec![start, start + 500]),
                Revalue::V(vec![start + 100, start + 900]),
            )
        } else {
            (Revalue::U(start), Revalue::U(start + 1000))
        };
        regions.push(Region {
            chrom: format!("chr{}", i % 5),
            start: s,
            end: e,
            score: String::from("."),
            strand: String::from(if i % 2 == 0 { "+" } else { "-" }),
            name: format!("r{}", i),
            regionlength: 1000,
        });
        let row: Vec<f32> = (0..ncols)
            .map(|_| match rng.next() % 10 {
                0 => f32::NAN,
                1 => 0.0,
                2 => -((rng.next() % 1000) as f32) / 7.0,
                _ => ((rng.next() % 100_000) as f32) / 3.0,
            })
            .collect();
        mat.push(row);
    }
    (regions, mat)
}

/// The previous serial row formatter (one String per value), as reference.
fn legacy_text(header: &str, regions: &[Region], mat: &[Vec<f32>], scale: f32) -> Vec<u8> {
    let mut t = header.to_string();
    for (region, row) in regions.iter().zip(mat.iter()) {
        t.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t",
            region.chrom,
            region.start.rewrites(),
            region.end.rewrites(),
            region.name,
            region.score,
            region.strand,
        ));
        t.push_str(
            &row.iter()
                .map(|x| {
                    if x.is_nan() {
                        "nan".to_string()
                    } else {
                        format!("{:.6}", scale * x)
                    }
                })
                .collect::<Vec<String>>()
                .join("\t"),
        );
        t.push('\n');
    }
    t.into_bytes()
}

#[test]
fn test_parallel_writer_matches_legacy_text() {
    let header = "@{\"upstream\":[3000],\"downstream\":[3000]}\n";
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(3)
        .build()
        .unwrap();
    for &(nrows, ncols, scale) in &[
        (0usize, 10usize, 1.0f32),
        (1, 0, 1.0),
        (1, 7, 1.0),
        (17, 3, 0.5),
        (30_000, 40, 1.0), // many chunks and several batches
        (500, 3000, 2.5),  // one row is larger than half a chunk
    ] {
        let (regions, mat) = test_regions_matrix(nrows, ncols, 42 + nrows as u64);
        let expected = legacy_text(header, &regions, &mat, scale);
        let mut gz: Vec<u8> = Vec::new();
        pool.install(|| {
            write_matrix_gz(
                &mut gz,
                header.as_bytes(),
                &mat,
                &regions,
                scale,
                Compression::default(),
            )
        })
        .unwrap();
        // GzDecoder (single member) also verifies CRC32 and ISIZE.
        let mut text = Vec::new();
        GzDecoder::new(&gz[..]).read_to_end(&mut text).unwrap();
        assert!(text == expected, "mismatch for {}x{}", nrows, ncols);
    }
}
