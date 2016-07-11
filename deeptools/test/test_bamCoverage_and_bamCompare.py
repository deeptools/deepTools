from nose.tools import assert_equal
import deeptools.bamCoverage as bam_cov
import deeptools.bamCompare as bam_comp
import deeptools.getScaleFactor as gs
import os.path
from os import unlink

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
BAMFILE_A = ROOT + "testA.bam"
BAMFILE_B = ROOT + "testB.bam"
BAMFILE_FILTER1 = ROOT + "test_filtering.bam"
BAMFILE_FILTER2 = ROOT + "test_filtering2.bam"
BEDFILE_FILTER = ROOT + "test_filtering.blacklist.bed"


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

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t0\t50\t0.00\n', '3R\t50\t150\t1.00\n', '3R\t150\t200\t2.00\n']
    assert_equal(resp, expected)
    unlink(outfile)


def test_bam_coverage_extend():
    outfile = '/tmp/test_file.bg'
    args = "-b {} -o {} --extendReads 100 --outFileFormat bedgraph".format(BAMFILE_B, outfile).split()
    bam_cov.main(args)
    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t0\t150\t1.00\n', '3R\t150\t200\t3.00\n']
    assert_equal(resp, expected)
    unlink(outfile)


def test_bam_coverage_extend_and_normalizeto1x():

    outfile = '/tmp/test_file.bg'
    args = "-b {} -o {} --normalizeTo1x 200 --extendReads 100 " \
           "--outFileFormat bedgraph".format(BAMFILE_B, outfile).split()
    bam_cov.main(args)
    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    # the scale factor should be 0.5, thus the result is similar to
    # that of the previous test divided by 0.5
    expected = ['3R\t0\t150\t0.50\n', '3R\t150\t200\t1.50\n']
    assert_equal(resp, expected)
    unlink(outfile)


def test_bam_coverage_skipnas():
    outfile = '/tmp/test_file.bg'
    args = "--bam {} -o {} --outFileFormat bedgraph --skipNAs".format(BAMFILE_B, outfile).split()
    bam_cov.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t50\t150\t1.00\n', '3R\t150\t200\t2.00\n']
    assert_equal(resp, expected)
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

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t0\t200\t1.00\n']
    assert_equal(resp, expected)
    unlink(outfile)


def test_bam_compare_diff_files():
    """
    Test with two different files
    """
    outfile = '/tmp/test_file.bg'
    args = "--bamfile1 {} --bamfile2 {} --scaleFactors 1:1 --ratio subtract " \
           "-o {} -p 1 --outFileFormat bedgraph".format(BAMFILE_A, BAMFILE_B, outfile).split()

    bam_comp.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t0\t50\t0.00\n', '3R\t50\t100\t-1.00\n', '3R\t100\t150\t0.00\n', '3R\t150\t200\t-1.00\n']
    assert_equal(resp, expected)
    unlink(outfile)


def test_get_num_kept_reads():
    """
    Test the scale factor functinos
    """
    args = "--bam {}  -o /tmp/test".format(BAMFILE_A, BAMFILE_B).split()

    args = bam_cov.process_args(args)
    num_kept_reads, total_reads = gs.get_num_kept_reads(args)

    # bam file 1 has 2 reads in 3R and 2 read in chr_cigar
    assert num_kept_reads == 3, "num_kept_reads is wrong"
    assert total_reads == 3, "num total reads is wrong"

    # ignore chr_cigar to count the total number of reads
    args = "--bam {} --ignoreForNormalization chr_cigar  -o /tmp/test".format(BAMFILE_A, BAMFILE_B).split()
    args = bam_cov.process_args(args)
    num_kept_reads, total_reads = gs.get_num_kept_reads(args)

    # the  number of kept reads should be 2 as the read on chr_cigar is skipped
    assert num_kept_reads == 2, "num_kept_reads is wrong ({})".format(num_kept_reads)

    # test filtering by read direction. Only forward reads are kept
    args = "--bam {}  -o /tmp/test --samFlagExclude 16 --ignoreForNormalization chr_cigar ".format(BAMFILE_A, BAMFILE_B).split()

    args = bam_cov.process_args(args)
    num_kept_reads, total_reads = gs.get_num_kept_reads(args)

    # only one forward read is expected in
    assert num_kept_reads == 1, "num_kept_reads is wrong"


def test_bam_compare_diff_files_skipnas():
    """
    Test skipnas
    Compared to the previous tests, any region that do not have coverage (in either of the bam files)
    is not included in the bedgraph file.
    """
    outfile = '/tmp/test_file.bg'
    args = "--bamfile1 {} --bamfile2 {} --scaleFactors 1:1 --ratio subtract " \
           "-o {} -p 1 --outFileFormat bedgraph --skipNAs".format(BAMFILE_A, BAMFILE_B, outfile).split()

    bam_comp.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t100\t150\t0.00\n', '3R\t150\t200\t-1.00\n']
    assert_equal(resp, expected)
    unlink(outfile)


def test_bam_compare_extend():
    """
    Test read extension
    """
    outfile = '/tmp/test_file.bg'
    args = "--bamfile1 {} --bamfile2 {} --extend 100 --scaleFactors 1:1 --ratio subtract " \
           "-o {} -p 1 --outFileFormat bedgraph".format(BAMFILE_A, BAMFILE_B, outfile).split()

    bam_comp.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t0\t100\t-1.00\n', '3R\t100\t150\t1.00\n', '3R\t150\t200\t-1.00\n']
    assert_equal(resp, expected)
    unlink(outfile)


def test_bam_compare_scale_factors_ratio():
    """
    Test scale factor
    """
    outfile = '/tmp/test_file.bg'
    args = "--bamfile1 {} --bamfile2 {} --ratio ratio --ignoreForNormalization chr_cigar " \
           "-o {} -p 1 --outFileFormat bedgraph".format(BAMFILE_A, BAMFILE_B, outfile).split()

    bam_comp.main(args)

    # The scale factors are [ 1.   0.5] because BAMFILE_B has dowble the amount of reads (4) compared to BAMFILE_A

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()

    """
    The distribution of reads for the bam file is:

                  0                              100                           200
                  |------------------------------------------------------------|
    testA.bam  3R                                ==============>
                                                                <==============


    testB.bam  3R                 <==============               ==============>
                                                 ==============>
                                                                ==============>

    ------------------------------------------------------------------------------

    ratio:             0      (0+1)/(1*0.5+1)=0.67             (1+1)/(1+2*0.5)=1
    (scale factors [1,0.5])                   (1+1)/(1+1*0.5)=1.33
    """

    expected = ['3R\t0\t50\t1.00\n', '3R\t50\t100\t0.67\n', '3R\t100\t150\t1.33\n', '3R\t150\t200\t1.00\n']
    assert_equal(resp, expected)
    unlink(outfile)


def test_bam_compare_scale_factors_subtract():
    """
    Test scale factor
    """
    outfile = '/tmp/test_file.bg'
    args = "--bamfile1 {} --bamfile2 {} --ratio subtract --ignoreForNormalization chr_cigar " \
           "-o {} -p 1 --outFileFormat bedgraph --normalizeTo1x 200".format(BAMFILE_A, BAMFILE_B, outfile).split()

    bam_comp.main(args)

    # The scale factors are [ 1.   0.5] because BAMFILE_B has dowble the amount of reads (4) compared to BAMFILE_A

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()

    """
    The distribution of reads for the bam file is:

                  0                              100                           200
                  |------------------------------------------------------------|
    testA.bam  3R                                ==============>
                                                                <==============


    testB.bam  3R                 <==============               ==============>
                                                 ==============>
                                                                ==============>

    ------------------------------------------------------------------------------

    subtract: scale factors [1,0.5], after applying normalize to 1x, coverage of test_A is 0.5, thus
    the factor to reach a coverate of 1 is x2. Thus, the final scale factors are [2,1]

    after applying factors:    0         -1              1              0

    """

    expected = ['3R\t0\t50\t0.00\n', '3R\t50\t100\t-1.00\n', '3R\t100\t150\t1.00\n', '3R\t150\t200\t0.00\n']
    assert_equal(resp, expected)
    unlink(outfile)


def test_bam_coverage_filter_blacklist():
    """
    Test --samFlagInclude --samFlagExclude --minMappingQuality --ignoreDuplicates and --blackListFileName
    """
    outfile = '/tmp/test_file_filter.bg'
    args = "--bam {} --normalizeTo1x 1400 -p 1 -o {} -of bedgraph --samFlagInclude 512 " \
           "--samFlagExclude 256 --minMappingQuality 5 --ignoreDuplicates " \
           "--blackListFileName {}".format(BAMFILE_FILTER1, outfile, BEDFILE_FILTER)
    args = args.split()
    bam_cov.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t0\t100\t0.00\n', '3R\t100\t150\t1.42\n', '3R\t150\t250\t4.88\n',
                '3R\t250\t300\t3.05\n', '3R\t300\t400\t2.24\n', '3R\t400\t450\t3.86\n',
                '3R\t450\t500\t4.07\n', '3R\t500\t550\t2.03\n', '3R\t550\t600\t2.44\n',
                '3R\t600\t650\t4.47\n', '3R\t650\t700\t3.46\n', '3R\t700\t750\t3.66\n',
                '3R\t750\t800\t4.07\n', '3R\t900\t950\t2.44\n', '3R\t950\t1000\t1.63\n',
                '3R\t1000\t1050\t0.81\n', '3R\t1050\t1500\t0.00\n']

    assert_equal(resp, expected)
    unlink(outfile)


def test_bam_compare_filter_blacklist():
    """
    Test --samFlagInclude --samFlagExclude --minMappingQuality --ignoreDuplicates and --blackListFileName
    """
    outfile = '/tmp/test_file_filter.bg'
    args = "-b1 {} -b2 {} --normalizeTo1x 1400 -p 1 -o {} -of bedgraph --samFlagInclude 512 " \
           "--samFlagExclude 256 --minMappingQuality 5 --ignoreDuplicates " \
           "--blackListFileName {}".format(BAMFILE_FILTER1, BAMFILE_FILTER2, outfile, BEDFILE_FILTER)
    args = args.split()
    bam_comp.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t0\t100\t0.00\n', '3R\t100\t150\t-0.22\n',
                '3R\t150\t200\t-0.16\n', '3R\t200\t250\t-0.07\n',
                '3R\t250\t300\t0.14\n', '3R\t300\t350\t0.10\n',
                '3R\t350\t400\t-0.09\n', '3R\t400\t450\t0.03\n',
                '3R\t450\t500\t0.10\n', '3R\t500\t550\t0.21\n',
                '3R\t550\t600\t0.02\n', '3R\t600\t650\t-0.10\n',
                '3R\t650\t700\t0.01\n', '3R\t700\t750\t-0.04\n',
                '3R\t750\t800\t-0.12\n', '3R\t900\t950\t0.21\n',
                '3R\t950\t1000\t0.20\n', '3R\t1000\t1050\t0.17\n',
                '3R\t1050\t1500\t0.00\n']
    assert_equal(resp, expected)
    unlink(outfile)
