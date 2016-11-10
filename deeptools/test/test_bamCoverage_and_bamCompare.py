from nose.tools import assert_equal
import deeptools.bamCoverage as bam_cov
import deeptools.bamCompare as bam_comp
import deeptools.getScaleFactor as gs
import os.path
import filecmp
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
    expected = ['3R\t0\t50\t0\n', '3R\t50\t150\t1\n', '3R\t150\t200\t2\n']
    assert_equal(resp, expected)
    unlink(outfile)


def test_bam_coverage_extend():
    outfile = '/tmp/test_file.bg'
    args = "-b {} -o {} --extendReads 100 --outFileFormat bedgraph".format(BAMFILE_B, outfile).split()
    bam_cov.main(args)
    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t0\t150\t1\n', '3R\t150\t200\t3\n']
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
    expected = ['3R\t0\t150\t0.5\n', '3R\t150\t200\t1.5\n']
    assert_equal(resp, expected)
    unlink(outfile)


def test_bam_coverage_skipnas():
    outfile = '/tmp/test_file.bg'
    args = "--bam {} -o {} --outFileFormat bedgraph --skipNAs".format(BAMFILE_B, outfile).split()
    bam_cov.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t50\t150\t1\n', '3R\t150\t200\t2\n']
    assert_equal(resp, expected)
    unlink(outfile)


def test_bam_coverage_filtering():
    outfile = '/tmp/test_file.bg'
    args = "--bam {} -o {} --outFileFormat bedgraph --ignoreDuplicates --verbose".format(BAMFILE_B, outfile).split()
    bam_cov.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t0\t50\t0\n', '3R\t50\t200\t1\n']
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
    expected = ['3R\t0\t200\t1\n']
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
    expected = ['3R\t0\t50\t0\n', '3R\t50\t100\t-1\n', '3R\t100\t150\t0\n', '3R\t150\t200\t-1\n']
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
    expected = ['3R\t100\t150\t0\n', '3R\t150\t200\t-1\n']
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
    expected = ['3R\t0\t100\t-1\n', '3R\t100\t150\t1\n', '3R\t150\t200\t-1\n']
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

    expected = ['3R\t0\t50\t1\n', '3R\t50\t100\t0.666667\n', '3R\t100\t150\t1.33333\n', '3R\t150\t200\t1\n']
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

    expected = ['3R\t0\t50\t0\n', '3R\t50\t100\t-1\n', '3R\t100\t150\t1\n', '3R\t150\t200\t0\n']
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
    expected = ['3R\t0\t100\t0\n', '3R\t100\t150\t1.42338\n', '3R\t150\t250\t4.88017\n',
                '3R\t250\t300\t3.05011\n', '3R\t300\t400\t2.23675\n', '3R\t400\t450\t3.86347\n',
                '3R\t450\t500\t4.06681\n', '3R\t500\t550\t2.03341\n', '3R\t550\t600\t2.44009\n',
                '3R\t600\t650\t4.47349\n', '3R\t650\t700\t3.45679\n', '3R\t700\t750\t3.66013\n',
                '3R\t750\t800\t4.06681\n', '3R\t900\t950\t2.44009\n', '3R\t950\t1000\t1.62672\n',
                '3R\t1000\t1050\t0.813362\n', '3R\t1050\t1500\t0\n']

    assert_equal(resp, expected)
    unlink(outfile)


def test_bam_coverage_offset1():
    """
    Test -bs 1 --Offset 1
    """
    outfile = '/tmp/test_offset.bw'
    args = "--Offset 1 --bam {} -p 1 -bs 1 -o {}".format(BAMFILE_A, outfile)
    args = args.split()
    bam_cov.main(args)
    try:
        # python 3 only
        filecmp.clear_cache()
    except:
        pass
    assert(filecmp.cmp(outfile, "{}testA_offset1.bw".format(ROOT)) is True)
    unlink(outfile)


def test_bam_coverage_offset1_10():
    """
    Test -bs 1 --Offset 1 10
    """
    outfile = '/tmp/test_offset.bw'
    args = "--Offset 1 10 -b {} -p 1 -bs 1 -o {}".format(BAMFILE_A, outfile)
    args = args.split()
    bam_cov.main(args)
    try:
        # python 3 only
        filecmp.clear_cache()
    except:
        pass
    assert(filecmp.cmp(outfile, "{}testA_offset1_10.bw".format(ROOT)) is True)
    unlink(outfile)


def test_bam_coverage_offset_minus1():
    """
    Test -bs 1 --Offset -1
    """
    outfile = '/tmp/test_offset.bw'
    args = "--Offset -1 -b {} -p 1 -bs 1 -o {}".format(BAMFILE_A, outfile)
    args = args.split()
    bam_cov.main(args)
    try:
        # python 3 only
        filecmp.clear_cache()
    except:
        pass
    assert(filecmp.cmp(outfile, "{}testA_offset-1.bw".format(ROOT)) is True)
    unlink(outfile)


def test_bam_coverage_offset20_minus4():
    """
    Test -bs 1 --Offset 20 -4
    """
    outfile = '/tmp/test_offset.bw'
    args = "--Offset 20 -4 -b {} -p 1 -bs 1 -o {}".format(BAMFILE_A, outfile)
    args = args.split()
    bam_cov.main(args)
    try:
        # python 3 only
        filecmp.clear_cache()
    except:
        pass
    assert(filecmp.cmp(outfile, "{}testA_offset20_-4.bw".format(ROOT)) is True)
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
    expected = ['3R\t0\t100\t0\n', '3R\t100\t150\t-0.220909\n',
                '3R\t150\t200\t-0.159356\n', '3R\t200\t250\t-0.0718929\n',
                '3R\t250\t300\t0.135883\n', '3R\t300\t350\t0.103093\n',
                '3R\t350\t400\t-0.0895516\n', '3R\t400\t450\t0.0308374\n',
                '3R\t450\t500\t0.0989418\n', '3R\t500\t550\t0.207044\n',
                '3R\t550\t600\t0.0198996\n', '3R\t600\t650\t-0.0957241\n',
                '3R\t650\t700\t0.00968255\n', '3R\t700\t750\t-0.040642\n',
                '3R\t750\t800\t-0.123451\n', '3R\t900\t950\t0.212545\n',
                '3R\t950\t1000\t0.199309\n', '3R\t1000\t1050\t0.167945\n',
                '3R\t1050\t1500\t0\n']
    assert_equal(resp, expected)
    unlink(outfile)
