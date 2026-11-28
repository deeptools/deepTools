use crate::covcalc::Revalue;
use crate::filehandler::{is_bed_or_gtf, read_bedfile, write_covfile};
use bigtools::Value;
use std::collections::HashMap;
use std::io::Write;
use tempfile::{NamedTempFile, tempdir};

fn write_bed(contents: &str) -> NamedTempFile {
    let mut f = NamedTempFile::new().expect("Failed to create temp BED file");
    write!(f, "{}", contents).expect("Failed to write temp BED file");
    f
}

fn chrom_chr1() -> HashMap<String, u32> {
    let mut chroms: HashMap<String, u32> = HashMap::new();
    chroms.insert("chr1".to_string(), 1_000_000);
    chroms
}

#[test]
fn test_read_bedfile_skips_start_ge_end_bed3() {
    let bed = write_bed("chr1\t100\t200\nchr1\t8000\t3000\nchr1\t300\t400\n");
    let chroms = chrom_chr1();

    let (regions, (_label, entries)) =
        read_bedfile(&bed.path().to_string_lossy().into_owned(), false, &chroms)
            .expect("read_bedfile failed");

    assert_eq!(entries, 2);
    assert_eq!(regions.len(), 2);
    assert_eq!(regions[0].start, Revalue::U(100));
    assert_eq!(regions[0].end, Revalue::U(200));
    assert_eq!(regions[0].regionlength, 100);
    assert_eq!(regions[1].start, Revalue::U(300));
    assert_eq!(regions[1].end, Revalue::U(400));
    assert_eq!(regions[1].regionlength, 100);
}

#[test]
fn test_read_bedfile_skips_start_ge_end_bed6() {
    let bed = write_bed(
        "chr1\t100\t200\tgeneA\t.\t+\nchr1\t13714\t11649\tgeneB\t.\t-\nchr1\t300\t400\tgeneC\t.\t-\n",
    );
    let chroms = chrom_chr1();

    let (regions, (_label, entries)) =
        read_bedfile(&bed.path().to_string_lossy().into_owned(), false, &chroms)
            .expect("read_bedfile failed");

    assert_eq!(entries, 2);
    assert_eq!(regions.len(), 2);
    assert_eq!(regions[0].name, "geneA");
    assert_eq!(regions[1].name, "geneC");
}

#[test]
fn test_read_bedfile_skips_start_eq_end() {
    let bed = write_bed("chr1\t100\t200\nchr1\t150\t150\n");
    let chroms = chrom_chr1();

    let (regions, (_label, entries)) =
        read_bedfile(&bed.path().to_string_lossy().into_owned(), false, &chroms)
            .expect("read_bedfile failed");

    assert_eq!(entries, 1);
    assert_eq!(regions.len(), 1);
}

#[test]
fn test_read_bedfile_missing_file_is_err() {
    let dir = tempdir().expect("Failed to create temp dir");
    let missing = dir.path().join("missing.bed");
    let chroms = chrom_chr1();

    let res = read_bedfile(&missing.to_string_lossy().into_owned(), false, &chroms);
    assert!(res.is_err());
}

#[test]
fn test_is_bed_or_gtf_missing_file_is_err() {
    let dir = tempdir().expect("Failed to create temp dir");
    let missing = dir.path().join("missing.bed");

    let res = is_bed_or_gtf(&missing.to_string_lossy().into_owned());
    assert!(res.is_err());
}

#[test]
fn test_read_bedfile_non_numeric_start_bed3_is_err() {
    let bed = write_bed("chr1\tabc\t200\n");
    let chroms = chrom_chr1();

    let err = read_bedfile(&bed.path().to_string_lossy().into_owned(), false, &chroms)
        .expect_err("non-numeric start should be an error");
    assert!(
        err.to_string().contains("chr1\tabc\t200"),
        "error should name the offending line, got: {err}"
    );
}

#[test]
fn test_read_bedfile_four_fields_is_err() {
    let bed = write_bed("chr1\t100\t200\tgeneA\n");
    let chroms = chrom_chr1();

    let res = read_bedfile(&bed.path().to_string_lossy().into_owned(), false, &chroms);
    assert!(res.is_err());
}

#[test]
fn test_write_covfile_bedgraph_missing_dir_is_err() {
    let dir = tempdir().expect("Failed to create temp dir");
    let ofile = dir.path().join("no_such_dir").join("out.bedgraph");
    let lines = std::iter::empty::<(String, Value)>();

    let res = write_covfile(
        lines,
        &ofile.to_string_lossy().into_owned(),
        "bedgraph",
        HashMap::new(),
    );
    assert!(res.is_err());
}
