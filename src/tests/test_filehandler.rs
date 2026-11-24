use crate::covcalc::Revalue;
use crate::filehandler::{read_bedfile, read_bedfile_skipped, unknown_chrom_warning};
use std::collections::HashMap;
use std::io::Write;
use tempfile::NamedTempFile;

fn write_bed(contents: &str) -> NamedTempFile {
    let mut f = NamedTempFile::new().expect("Failed to create temp BED file");
    write!(f, "{}", contents).expect("Failed to write temp BED file");
    f
}

#[test]
fn test_read_bedfile_skips_start_ge_end_bed3() {
    let bed = write_bed("chr1\t100\t200\nchr1\t8000\t3000\nchr1\t300\t400\n");
    let mut chroms: HashMap<String, u32> = HashMap::new();
    chroms.insert("chr1".to_string(), 1_000_000);

    let (regions, (_label, entries)) =
        read_bedfile(&bed.path().to_string_lossy().into_owned(), false, &chroms);

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
    let mut chroms: HashMap<String, u32> = HashMap::new();
    chroms.insert("chr1".to_string(), 1_000_000);

    let (regions, (_label, entries)) =
        read_bedfile(&bed.path().to_string_lossy().into_owned(), false, &chroms);

    assert_eq!(entries, 2);
    assert_eq!(regions.len(), 2);
    assert_eq!(regions[0].name, "geneA");
    assert_eq!(regions[1].name, "geneC");
}

#[test]
fn test_read_bedfile_skips_start_eq_end() {
    let bed = write_bed("chr1\t100\t200\nchr1\t150\t150\n");
    let mut chroms: HashMap<String, u32> = HashMap::new();
    chroms.insert("chr1".to_string(), 1_000_000);

    let (regions, (_label, entries)) =
        read_bedfile(&bed.path().to_string_lossy().into_owned(), false, &chroms);

    assert_eq!(entries, 1);
    assert_eq!(regions.len(), 1);
}

#[test]
fn test_read_bedfile_counts_unknown_chroms() {
    let bed = write_bed("chr1\t100\t200\nchrZ\t1\t50\nchrZ\t60\t90\nchrY\t5\t10\n");
    let mut chroms: HashMap<String, u32> = HashMap::new();
    chroms.insert("chr1".to_string(), 1_000_000);

    let (regions, (_label, entries), unknown) =
        read_bedfile_skipped(&bed.path().to_string_lossy().into_owned(), false, &chroms);

    assert_eq!(entries, 1);
    assert_eq!(regions.len(), 1);
    assert_eq!(unknown.get("chrZ"), Some(&2));
    assert_eq!(unknown.get("chrY"), Some(&1));
    assert_eq!(unknown.len(), 2);
}

#[test]
fn test_read_bedfile_empty_chroms_skips_everything_and_says_why() {
    let bed = write_bed("chr1\t100\t200\nchr2\t1\t50\n");
    let chroms: HashMap<String, u32> = HashMap::new();
    let path = bed.path().to_string_lossy().into_owned();

    let (regions, (_label, entries), unknown) = read_bedfile_skipped(&path, false, &chroms);

    assert!(regions.is_empty());
    assert_eq!(entries, 0);
    let msg = unknown_chrom_warning(&path, &chroms, &unknown).expect("expected a warning");
    assert!(msg.contains("no chromosome sizes"), "got: {msg}");
}

#[test]
fn test_unknown_chrom_warning_names_chroms_and_counts() {
    let bed = write_bed("chr1\t100\t200\nchrZ\t1\t50\nchrZ\t60\t90\n");
    let mut chroms: HashMap<String, u32> = HashMap::new();
    chroms.insert("chr1".to_string(), 1_000_000);
    let path = bed.path().to_string_lossy().into_owned();

    let (_, _, unknown) = read_bedfile_skipped(&path, false, &chroms);
    let msg = unknown_chrom_warning(&path, &chroms, &unknown).expect("expected a warning");

    assert!(msg.contains("chrZ (2)"), "got: {msg}");
    assert!(!msg.contains("chr1"), "got: {msg}");
}

#[test]
fn test_unknown_chrom_warning_silent_when_all_known() {
    let bed = write_bed("chr1\t100\t200\n");
    let mut chroms: HashMap<String, u32> = HashMap::new();
    chroms.insert("chr1".to_string(), 1_000_000);
    let path = bed.path().to_string_lossy().into_owned();

    let (_, _, unknown) = read_bedfile_skipped(&path, false, &chroms);

    assert!(unknown_chrom_warning(&path, &chroms, &unknown).is_none());
}
