use crate::covcalc::Revalue;
use crate::filehandler::{read_bedfile, write_covfile};
use bigtools::{BigWigRead, Value};
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

fn chrom_ordered_lines(order: &[&str]) -> Vec<(String, Value)> {
    order
        .iter()
        .map(|c| {
            (
                c.to_string(),
                Value {
                    start: 0,
                    end: 100,
                    value: 1.0,
                },
            )
        })
        .collect()
}

#[test]
fn test_write_covfile_bigwig_with_non_lexicographic_chrom_order() {
    // BAM headers list chromosomes as chr1, chr2, chr10; bigWig input must be sorted by name.
    let sizes: HashMap<String, u32> = [("chr1", 1000), ("chr2", 1000), ("chr10", 1000)]
        .iter()
        .map(|(c, l)| (c.to_string(), *l))
        .collect();
    let out = NamedTempFile::new().expect("temp bigwig");
    let path = out.path().to_string_lossy().into_owned();

    write_covfile(
        chrom_ordered_lines(&["chr1", "chr2", "chr10"]).into_iter(),
        &path,
        "bigwig",
        sizes,
    );

    let mut bw = BigWigRead::open_file(&path).expect("bigwig must be readable");
    for chrom in ["chr1", "chr2", "chr10"] {
        let intervals: Vec<Value> = bw
            .get_interval(chrom, 0, 1000)
            .unwrap_or_else(|e| panic!("{chrom}: {e}"))
            .collect::<Result<_, _>>()
            .unwrap_or_else(|e| panic!("{chrom}: unreadable interval block: {e}"));
        assert_eq!(intervals.len(), 1, "{chrom} lost its interval");
        assert_eq!(intervals[0].start, 0, "{chrom} start");
        assert_eq!(intervals[0].end, 100, "{chrom} end");
        assert_eq!(intervals[0].value, 1.0, "{chrom} value");
    }
}
