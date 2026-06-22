import deeptools.alignmentSieve as aln_seive

import os.path
from os import unlink

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
BAMFILE_IN = ROOT + "paired_chr2L.bam"
CRAMFILE_IN = ROOT + "paired_chr2L.cram"
FILTER_METRICS_FILE = ROOT + "alignmentSieve.txt"
BED_OUT = ROOT + "alignmentSieve.bed"
BAMFILE_OUT = ROOT + "alignmentSieve.bam"
BAMFILE_OUT2 = ROOT + "alignmentSieve2.bam"
BAMFILE_OUT3 = ROOT + "alignmentSieve3.bam"


def test_minimum_mapping_quality_filter_metric():
    """
    Test minimal command line args for alignement sieve
    """
    outfile = '/tmp/test_bam.bam'
    filter_metric_file = '/tmp/test_metrics.txt'
    args = "--bam {} -o {} --minMappingQuality 10 --filterMetrics {}".format(BAMFILE_IN, outfile, filter_metric_file).split()
    aln_seive.main(args)

    _foo = open(filter_metric_file, "r")
    resp = _foo.readlines()[2]
    _foo.close()
    expected = 'paired_chr2L.bam\t8440\t12644\n'
    assert expected in resp, f"'{expected}' not found in '{resp}'"

    bam_file_size = os.path.getsize(BAMFILE_OUT)
    expected_file_size = os.path.getsize(outfile)
    size_tolerance = 5000
    size_difference = abs(bam_file_size - expected_file_size)
    assert size_difference <= size_tolerance, "File size do not match"
    unlink(outfile)
    unlink(filter_metric_file)


def test_with_bed_output_along_with_shifts():
    """
    Tests Alignment seive with shifts and output BED file
    """
    output_bed_file = "/tmp/aln_seive.bed"
    args = "--bam {} -o {} --minMappingQuality 10 --BED --shift 1 -2 3 -4".format(BAMFILE_IN, output_bed_file).split()
    aln_seive.main(args)
    with open(output_bed_file, "r") as _foo:
        result = len(_foo.readlines())
    _expected = 4261
    assert result == _expected, "No of lines in BED files differ"
    unlink(output_bed_file)


def test_with_bam_output_with_shifts():
    """
    Tests Alignment seive with shifts and output BAM file
    """
    output_bam_file = "/tmp/aln_seive.bam"
    args = "--bam {} -o {} --minMappingQuality 10 --shift 1 -2 3 -4".format(BAMFILE_IN, output_bam_file).split()
    aln_seive.main(args)

    bam_file_size = os.path.getsize(BAMFILE_OUT2)
    expected_file_size = os.path.getsize(output_bam_file)
    size_tolerance = 5000
    size_difference = abs(bam_file_size - expected_file_size)
    assert size_difference <= size_tolerance, "File sizes do not match"
    unlink(output_bam_file)


def test_with_cram_output_with_shifts():
    """
    Tests Alignment seive with  CRAM input along with shifts
    """
    output_bam_file = "/tmp/aln_seive2.bam"
    args = "--bam {} -o {} --minMappingQuality 10 --shift 1 -2 3 -4".format(BAMFILE_IN, output_bam_file).split()
    aln_seive.main(args)

    bam_file_size = os.path.getsize(BAMFILE_OUT3)
    expected_file_size = os.path.getsize(output_bam_file)
    size_tolerance = 5000
    size_difference = abs(bam_file_size - expected_file_size)
    assert size_difference <= size_tolerance, "File sizes do not match"
    unlink(output_bam_file)
