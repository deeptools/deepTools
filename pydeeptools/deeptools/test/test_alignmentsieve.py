import deeptools.alignmentSieve2 as aln_seive
import tempfile
from pathlib import Path
import pysam

ROOT = Path(__file__).parent / "test_data"
BAMFILE_IN = str(ROOT / "paired_chr2L.bam")
CRAMFILE_IN = str(ROOT / "paired_chr2L.cram")
BL_IN = str(ROOT / "alignmentSieve_blacklist.bed")
FILTER_METRICS_FILE = str(ROOT / "alignmentSieve.txt")
BED_OUT = str(ROOT / "alignmentSieve.bed")
BEDFILE = str(ROOT / "alignmentSieve.bed")
BAMFILE_OUT1 = str(ROOT / "alignmentSieve1.bam")
BAMFILEFILT_OUT1 = str(ROOT / "alignmentSieve1.filtered.bam")
BAMFILE_OUT2 = str(ROOT / "alignmentSieve2.bam")
BAMFILE_OUT3 = str(ROOT / "alignmentSieve3.bam")
BAMFILE_OUT4 = str(ROOT / "alignmentSieve4.bam")
BAMFILE_OUT5 = str(ROOT / "alignmentSieve5.bam")
BAMFILE_OUT6 = str(ROOT / "alignmentSieve6.bam")
BAMFILE_OUT7 = str(ROOT / "alignmentSieve7.bam")

def _assert_equals(expected: str, actual: str) -> None:
    if Path(expected).suffix == '.bam':
        with (
            pysam.AlignmentFile(expected, "rb") as ef,
            pysam.AlignmentFile(actual, "rb") as af,
        ):
            for eread, aread in zip(ef, af):
                assert eread.to_string() == aread.to_string(), f"BAM files do not match: {eread} vs {aread}"
    else:
        with (
            open(expected, 'r') as ef,
            open(actual, 'r') as af,
        ):
            for eread, aread in zip(ef, af):
                assert eread == aread, f"files do not match: {eread} vs {aread}"

def test_alsieve_minmapq():
    """
    Test minimal command line args for alignement sieve
    """
    _, outfile = tempfile.mkstemp(suffix=".bam")
    _, filter_metric_file = tempfile.mkstemp(suffix=".txt")
    _, filtered_bamfile = tempfile.mkstemp(suffix=".bam")
    args = f"--bam {BAMFILE_IN} -o {outfile} --minMappingQuality 10 --filteredOutReads {filtered_bamfile} --filterMetrics {filter_metric_file}".split()
    aln_seive.main(args)

    _foo = open(filter_metric_file, "r")
    resp = _foo.readlines()[2]
    _foo.close()
    expected = 'paired_chr2L.bam\t8440\t12644\n'
    assert expected in resp, f"'{expected}' not found in '{resp}'"

    _assert_equals(BAMFILE_OUT1, outfile)


def test_alsieve_fraglen():
    """
    Test --minFragmentLength / --maxFragmentLength filtering (fragment-size
    selection, e.g. nucleosome-free vs mono-nucleosome). Only read pairs whose
    template length falls in [130, 200] are kept.
    """
    _, outfile = tempfile.mkstemp(suffix=".bam")
    _, filter_metric_file = tempfile.mkstemp(suffix=".txt")
    args = f"--bam {BAMFILE_IN} -o {outfile} --minFragmentLength 130 --maxFragmentLength 200 --filterMetrics {filter_metric_file}".split()
    aln_seive.main(args)

    _foo = open(filter_metric_file, "r")
    resp = _foo.readlines()[2]
    _foo.close()
    expected = 'paired_chr2L.bam\t7933\t12644\n'
    assert expected in resp, f"'{expected}' not found in '{resp}'"
    _assert_equals(BAMFILE_OUT2, outfile)


def test_alsieve_rnastrand():
    """
    Test --filterRNAstrand forward (keep only reads matching the forward
    transcription strand, using paired-end mate orientation).
    """
    _, outfile = tempfile.mkstemp(suffix=".bam")
    _, filter_metric_file = tempfile.mkstemp(suffix=".txt")
    args = f"--bam {BAMFILE_IN} -o {outfile} --filterRNAstrand forward --filterMetrics {filter_metric_file}".split()
    aln_seive.main(args)

    _foo = open(filter_metric_file, "r")
    resp = _foo.readlines()[2]
    _foo.close()
    expected = 'paired_chr2L.bam\t6303\t12644\n'
    assert expected in resp, f"'{expected}' not found in '{resp}'"
    _assert_equals(BAMFILE_OUT3, outfile)


def test_alsieve_shift_bed():
    """
    Tests Alignment seive with shifts and output BED file
    """
    _, output_bed_file = tempfile.mkstemp(suffix=".bed")
    args = f"--bam {BAMFILE_IN} -o {output_bed_file} --minMappingQuality 10 --BED --shift 1 -2 3 -4".split()
    aln_seive.main(args)
    with open(output_bed_file, "r") as _foo:
        result = len(_foo.readlines())
    _expected = 4261
    assert result == _expected, "No of lines in BED files differ"
    _assert_equals(BEDFILE, output_bed_file)


def test_alsieve_shift():
    """
    Tests Alignment seive with shifts and output BAM file
    """
    _, output_bam_file = tempfile.mkstemp(suffix=".bam")
    args = f"--bam {BAMFILE_IN} -o {output_bam_file} --minMappingQuality 10 --shift 1 -2 3 -4".split()
    aln_seive.main(args)

    _assert_equals(BAMFILE_OUT4, output_bam_file)


def test_alsieve_blacklist():
    """
    Tests Alignment seive with shifts and output BAM file
    """
    _, output_bam_file = tempfile.mkstemp(suffix=".bam")
    args = f"--bam {BAMFILE_IN} -bl {BL_IN} -o {output_bam_file} --minMappingQuality 10 --shift 1 -2 3 -4".split()
    aln_seive.main(args)

    _assert_equals(BAMFILE_OUT5, output_bam_file)

def test_alsieve_shift2():
    """
    shift with 2 arguments.
    """
    _, output_bam_file = tempfile.mkstemp(suffix=".bam")
    args = f"--bam {BAMFILE_IN} -o {output_bam_file} --shift 1 -2".split()
    aln_seive.main(args)

    _assert_equals(BAMFILE_OUT6, output_bam_file)

def test_alsieve_atacshift():
    """
    Tests Alignment seive with ATAC shift
    """
    _, output_bam_file = tempfile.mkstemp(suffix=".bam")
    args = f"--bam {BAMFILE_IN} -o {output_bam_file} --minMappingQuality 10 --ATACshift".split()
    aln_seive.main(args)

    _assert_equals(BAMFILE_OUT7, output_bam_file)
