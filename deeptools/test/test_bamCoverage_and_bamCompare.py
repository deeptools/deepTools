import deeptools.bamCoverage as bam_cov
import deeptools.bamCompare as bam_comp
import os.path
from os import unlink

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
BAMFILE_A = ROOT + "testA.bam"
BAMFILE_B = ROOT + "testB.bam"


"""
The distribution of reads for the bam file is:

              0                              100                           200
              |------------------------------------------------------------|
testA.bam  3R                                ==============>
                                                            <==============


testB.bam  3R                 <==============               ==============>
                                             ==============>
                                                            ==============>
        """


def test_bam_coverage_arguments():
    """
    Test minimal command line args for bamCoverage
    """
    outfile = '/tmp/test_file.bg'
    args = "--bam {} -o {} --outFileFormat bedgraph".format(BAMFILE_B, outfile).split()
    bam_cov.main(args)

    resp = open(outfile, 'r').readlines()
    expected = ['3R\t0\t50\t0.00\n', '3R\t50\t150\t1.00\n', '3R\t150\t200\t2.0\n']
    assert resp == expected, "{} != {}".format(resp, expected)
    unlink(outfile)


def test_bam_coverage_extend():
    outfile = '/tmp/test_file.bg'
    args = "-b {} -o {} --extendReads 100 --outFileFormat bedgraph".format(BAMFILE_B, outfile).split()
    bam_cov.main(args)
    resp = open(outfile, 'r').readlines()

    assert resp == ['3R\t0\t150\t1.00\n', '3R\t150\t200\t3.0\n']
    unlink(outfile)


def test_bam_coverage_extend_and_normalizeto1x():

    outfile = '/tmp/test_file.bg'
    args = "-b {} -o {} --normalizeTo1x 200 --extendReads 100 " \
           "--outFileFormat bedgraph".format(BAMFILE_B, outfile).split()
    bam_cov.main(args)
    resp = open(outfile, 'r').readlines()
    # the scale factor should be 0.5, thus the result is similar to
    # that of the previous test divided by 0.5
    expected = ['3R\t0\t150\t0.50\n', '3R\t150\t200\t1.5\n']
    assert resp == expected, "{} != {}".format(resp, expected)
    unlink(outfile)


def test_bam_coverage_skipnas():
    outfile = '/tmp/test_file.bg'
    args = "--bam {} -o {} --outFileFormat bedgraph --skipNAs".format(BAMFILE_B, outfile).split()
    bam_cov.main(args)

    resp = open(outfile, 'r').readlines()
    expected = ['3R\t50\t150\t1.00\n', '3R\t150\t200\t2.0\n']
    assert resp == expected, "{} != {}".format(resp, expected)
    unlink(outfile)


def test_bam_compare_arguments():
    """
    Test minimal command line args for bamCoverage. The ratio
    between the same file is taken, therefore, the expected value
    is 1.0 for all bins.
    """
    outfile = '/tmp/test_file.bg'
    args = "--bamfile1 {} --bamfile2 {} " \
           "-o {} -p 1 --outFileFormat bedgraph --ratio ratio".format(BAMFILE_B, BAMFILE_B, outfile).split()

    bam_comp.main(args)

    resp = open(outfile, 'r').readlines()
    expected = ['3R\t0\t200\t1.0\n']
    assert resp == expected, "{} != {}".format(resp, expected)
    unlink(outfile)


def test_bam_compare_diff_files():
    """
    Test with two different files
    """
    outfile = '/tmp/test_file.bg'
    args = "--bamfile1 {} --bamfile2 {} --scaleFactors 1:1 --ratio subtract " \
           "-o {} -p 1 --outFileFormat bedgraph".format(BAMFILE_A, BAMFILE_B, outfile).split()

    bam_comp.main(args)

    resp = open(outfile, 'r').readlines()
    expected = ['3R\t0\t50\t0.00\n', '3R\t50\t100\t-1.00\n', '3R\t100\t150\t0.00\n', '3R\t150\t200\t-1.0\n']
    assert resp == expected, "{} != {}".format(resp, expected)
    unlink(outfile)


def test_bam_compare_diff_files_skipnas():
    """
    Test skipnas
    Compared to the previous tests, any region that do not have coverage (in either of the bam files)
    is not included in the bedgraphfile.
    """
    outfile = '/tmp/test_file.bg'
    args = "--bamfile1 {} --bamfile2 {} --scaleFactors 1:1 --ratio subtract " \
           "-o {} -p 1 --outFileFormat bedgraph --skipNAs".format(BAMFILE_A, BAMFILE_B, outfile).split()

    bam_comp.main(args)

    resp = open(outfile, 'r').readlines()
    expected = ['3R\t100\t150\t0.00\n', '3R\t150\t200\t-1.0\n']
    assert resp == expected, "{} != {}".format(resp, expected)
    unlink(outfile)


def test_bam_compare_extend():
    """
    Test read extension
    """
    outfile = '/tmp/test_file.bg'
    args = "--bamfile1 {} --bamfile2 {} --extend 100 --scaleFactors 1:1 --ratio subtract " \
           "-o {} -p 1 --outFileFormat bedgraph".format(BAMFILE_A, BAMFILE_B, outfile).split()

    bam_comp.main(args)

    resp = open(outfile, 'r').readlines()
    expected = ['3R\t0\t100\t-1.00\n', '3R\t100\t150\t1.00\n', '3R\t150\t200\t-1.0\n']
    assert resp == expected, "{} != {}".format(resp, expected)
    unlink(outfile)
