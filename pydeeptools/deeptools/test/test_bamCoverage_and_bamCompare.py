import deeptools.bamCoverage2 as bam_cov
import deeptools.bamCompare2 as bam_comp
import os.path
import filecmp
from os import unlink

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
BAMFILE_A = ROOT + "testA.bam"
BAMFILE_B = ROOT + "testB.bam"
BAMFILE_FILTER1 = ROOT + "test_filtering.bam"
BAMFILE_FILTER2 = ROOT + "test_filtering2.bam"
CRAMFILE_A = ROOT + "testA.cram"
CRAMFILE_B = ROOT + "testB.cram"
CRAMFILE_FILTER1 = ROOT + "test_filtering.cram"
CRAMFILE_FILTER2 = ROOT + "test_filtering2.cram"
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
    #for fname in [BAMFILE_B, CRAMFILE_B]:
    for fname in [BAMFILE_B]:
        args = "--bam {} -o {} --outFileFormat bedgraph".format(fname, outfile).split()
        bam_cov.main(args)

        _foo = open(outfile, 'r')
        resp = _foo.readlines()
        _foo.close()
        expected = ['3R\t0\t50\t0\n', '3R\t50\t150\t1\n', '3R\t150\t200\t2\n']
        assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
        unlink(outfile)


def test_bam_coverage_extend():
    outfile = '/tmp/test_file.bg'
    #for fname in [BAMFILE_B, CRAMFILE_B]:
    for fname in [BAMFILE_B]:
        args = "-b {} -o {} --extendReads 100 --outFileFormat bedgraph".format(fname, outfile).split()
        bam_cov.main(args)
        _foo = open(outfile, 'r')
        resp = _foo.readlines()
        _foo.close()
        expected = ['3R\t0\t150\t1\n', '3R\t150\t200\t3\n']
        assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
        unlink(outfile)


def test_bam_coverage_extend_and_normalizeUsingRPGC():

    outfile = '/tmp/test_file.bg'
    #for fname in [BAMFILE_B, CRAMFILE_B]:
    for fname in [BAMFILE_B]:
        args = "-b {} -o {} --normalizeUsing RPGC --effectiveGenomeSize 200 --extendReads 100 --verbose " \
               "--outFileFormat bedgraph".format(fname, outfile).split()
        bam_cov.main(args)
        _foo = open(outfile, 'r')
        resp = _foo.readlines()
        _foo.close()
        # the scale factor should be 0.5, thus the result is similar to
        # that of the previous test divided by 0.5
        expected = ['3R\t0\t150\t0.5\n', '3R\t150\t200\t1.5\n']
        assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
        unlink(outfile)


def test_bam_coverage_skipnas():
    outfile = '/tmp/test_file.bg'
    #for fname in [BAMFILE_B, CRAMFILE_B]:
    for fname in [BAMFILE_B]:
        args = "--bam {} -o {} --outFileFormat bedgraph --skipNAs".format(fname, outfile).split()
        bam_cov.main(args)

        _foo = open(outfile, 'r')
        resp = _foo.readlines()
        _foo.close()
        expected = ['3R\t50\t150\t1\n', '3R\t150\t200\t2\n']
        assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
        unlink(outfile)


# def test_bam_coverage_filtering():
#     outfile = '/tmp/test_file.bg'
#     #for fname in [BAMFILE_B, CRAMFILE_B]:
#     for fname in [BAMFILE_B]:
#         args = "--bam {} -o {} --outFileFormat bedgraph --ignoreDuplicates --verbose".format(fname, outfile).split()
#         bam_cov.main(args)

#         _foo = open(outfile, 'r')
#         resp = _foo.readlines()
#         _foo.close()
#         expected = ['3R\t0\t50\t0\n', '3R\t50\t200\t1\n']
#         assert resp == expected, "{} != {}".format(resp, expected)
#         unlink(outfile)


def test_bam_compare_arguments():
    """
    Test minimal command line args for bamCoverage. The ratio
    between the same file is taken, therefore, the expected value
    is 1.0 for all bins.
    """
    outfile = '/tmp/test_file.bg'
    #for fname in [BAMFILE_B, CRAMFILE_B]:
    for fname in [BAMFILE_B]:
        args = "--bamfile1 {} --bamfile2 {} " \
               "-o {} -p 1 --outFileFormat bedgraph --operation ratio".format(fname, fname, outfile).split()
        bam_comp.main(args)

        _foo = open(outfile, 'r')
        resp = _foo.readlines()
        _foo.close()
        expected = ['3R\t0\t200\t1\n']
        assert resp == expected, "{} != {}".format(resp, expected)
        unlink(outfile)


def test_bam_compare_diff_files():
    """
    Test with two different files
    """
    outfile = '/tmp/test_file.bg'
    #for A, B in [(BAMFILE_A, BAMFILE_B), (CRAMFILE_A, CRAMFILE_B)]:
    for A, B in [(BAMFILE_A, BAMFILE_B)]:
        args = "--bamfile1 {} --bamfile2 {} --scaleFactors 1:1 --operation subtract --verbose " \
               "-o {} -p 1 --outFileFormat bedgraph".format(A, B, outfile).split()
        bam_comp.main(args)

        _foo = open(outfile, 'r')
        resp = _foo.readlines()
        _foo.close()
        expected = ['3R\t0\t50\t0\n', '3R\t50\t100\t-1\n', '3R\t100\t150\t0\n', '3R\t150\t200\t-1\n']
        assert resp == expected, "{} != {}".format(resp, expected)
        unlink(outfile)


def test_bam_compare_pseudocounts():
    """
    Test with different pseudocounts
    """
    outfile = '/tmp/test_file.bg'
    args = "--bamfile1 {} --bamfile2 {} --outFileFormat bedgraph --scaleFactors 1:1 -o {} " \
           "--pseudocount 1 0".format(BAMFILE_A, BAMFILE_B, outfile).split()
    bam_comp.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t0\t50\tinf\n', '3R\t50\t100\t0\n', '3R\t100\t150\t1\n', '3R\t150\t200\t0\n']
    assert resp == expected, "{} != {}".format(resp, expected)
    unlink(outfile)


def test_bam_compare_ZoverZ():
    """
    Ensure --skipZeroOverZero works in bamCompare
    """
    outfile = '/tmp/test_file.bg'
    args = "--bamfile1 {} --bamfile2 {} --outFileFormat bedgraph --scaleFactors 1:1 -o {} " \
           "--skipZeroOverZero --verbose".format(BAMFILE_A, BAMFILE_B, outfile).split()
    bam_comp.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t50\t100\t-1\n', '3R\t100\t150\t0\n', '3R\t150\t200\t-0.58\n']
    assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
    unlink(outfile)


# def test_get_num_kept_reads():
#     """
#     Test the scale factor functions
#     """
#     for fname in [BAMFILE_A, CRAMFILE_A]:
#         args = "--bam {}  -o /tmp/test".format(fname).split()

#         args = bam_cov.process_args(args)
#         num_kept_reads, total_reads = gs.get_num_kept_reads(args, None)

#         # bam file 1 has 2 reads in 3R and 2 read in chr_cigar
#         assert num_kept_reads == 3, "num_kept_reads is wrong"
#         assert total_reads == 3, "num total reads is wrong"

#         # ignore chr_cigar to count the total number of reads
#         args = "--bam {} --ignoreForNormalization chr_cigar  -o /tmp/test".format(fname).split()
#         args = bam_cov.process_args(args)
#         num_kept_reads, total_reads = gs.get_num_kept_reads(args, None)

#         # the  number of kept reads should be 2 as the read on chr_cigar is skipped
#         assert num_kept_reads == 2, "num_kept_reads is wrong ({})".format(num_kept_reads)

#         # test filtering by read direction. Only forward reads are kept
#         args = "--bam {}  -o /tmp/test --samFlagExclude 16 --ignoreForNormalization chr_cigar ".format(fname).split()

#         args = bam_cov.process_args(args)
#         num_kept_reads, total_reads = gs.get_num_kept_reads(args, None)

#         # only one forward read is expected in
#         assert num_kept_reads == 1, "num_kept_reads is wrong"


def test_bam_compare_diff_files_skipnas():
    """
    Test skipnas
    Compared to the previous tests, any region that do not have coverage (in either of the bam files)
    is not included in the bedgraph file.
    """
    outfile = '/tmp/test_file.bg'
    #for A, B in [(BAMFILE_A, BAMFILE_B), (CRAMFILE_A, CRAMFILE_B)]:
    for A, B in [(BAMFILE_A, BAMFILE_B)]:
        args = "--bamfile1 {} --bamfile2 {} --scaleFactors 1:1 --operation subtract " \
               "-o {} -p 1 --outFileFormat bedgraph --skipNAs".format(A, B, outfile).split()
        bam_comp.main(args)

        _foo = open(outfile, 'r')
        resp = _foo.readlines()
        _foo.close()
        expected = ['3R\t100\t150\t0\n', '3R\t150\t200\t-1\n']
        assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
        unlink(outfile)


def test_bam_compare_extend():
    """
    Test read extension
    """
    outfile = '/tmp/test_file.bg'
    #for A, B in [(BAMFILE_A, BAMFILE_B), (CRAMFILE_A, CRAMFILE_B)]:
    for A, B in [(BAMFILE_A, BAMFILE_B)]:
        args = "--bamfile1 {} --bamfile2 {} --extend 100 --scaleFactors 1:1 --operation subtract " \
               "-o {} -p 1 --outFileFormat bedgraph".format(A, B, outfile).split()
        bam_comp.main(args)

        _foo = open(outfile, 'r')
        resp = _foo.readlines()
        _foo.close()
        expected = ['3R\t0\t100\t-1\n', '3R\t100\t150\t1\n', '3R\t150\t200\t-1\n']
        assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
        unlink(outfile)


def test_bam_compare_scale_factors_ratio():
    """
    Test scale factor
    """
    outfile = '/tmp/test_file.bg'
    #for A, B in [(BAMFILE_A, BAMFILE_B), (CRAMFILE_A, CRAMFILE_B)]:
    for A, B in [(BAMFILE_A, BAMFILE_B)]:
        args = "--bamfile1 {} --bamfile2 {} --operation ratio --ignoreForNormalization chr_cigar " \
               "-o {} -p 1 --outFileFormat bedgraph".format(A, B, outfile).split()
        bam_comp.main(args)

        # The scale factors are [ 1.   0.5] because BAMFILE_B has double the amount of reads (4) compared to BAMFILE_A

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

        expected = ['3R\t0\t50\t1\n', '3R\t50\t100\t0.67\n', '3R\t100\t150\t1.33\n', '3R\t150\t200\t1\n']
        assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
        unlink(outfile)


def test_bam_compare_scale_factors_subtract():
    """
    Test scale factor
    """
    outfile = '/tmp/test_file.bg'
    #for A, B in [(BAMFILE_A, BAMFILE_B), (CRAMFILE_A, CRAMFILE_B)]:
    for A, B in [(BAMFILE_A, BAMFILE_B)]:
        args = "--bamfile1 {} --bamfile2 {} --operation subtract --ignoreForNormalization chr_cigar " \
               "-o {} -p 1 --outFileFormat bedgraph --scaleFactorsMethod None --normalizeUsing CPM --verbose".format(A, B, outfile).split()
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

        subtract: After applying CPM normalization, the scale factors are [500000,250000]

        after applying factors:    0         -25k              25k              0

        """

        expected = ['3R\t0\t50\t0\n', '3R\t50\t100\t-250000\n', '3R\t100\t150\t250000\n', '3R\t150\t200\t0\n']
        assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
        unlink(outfile)


def test_bam_coverage_filter_blacklist():
    """
    Test --samFlagInclude --samFlagExclude --minMappingQuality and --blackListFileName
    """
    outfile = '/tmp/test_file_filter.bg'
    #for fname in [BAMFILE_FILTER1, CRAMFILE_FILTER1]:
    for fname in [BAMFILE_FILTER1]:
        args = "--bam {} --normalizeUsing RPGC --effectiveGenomeSize 1400 -p 1 -o {} -of bedgraph --samFlagInclude 512 " \
               "--samFlagExclude 256 --minMappingQuality 5 --verbose " \
               "--blackListFileName {}".format(fname, outfile, BEDFILE_FILTER)
        args = args.split()
        bam_cov.main(args)

        _foo = open(outfile, 'r')
        resp = _foo.readlines()
        _foo.close()

        expected = [
            "3R\t0\t100\t0\n",
            "3R\t100\t150\t1.75\n",
            "3R\t150\t200\t5.84\n",
            "3R\t200\t250\t6.04\n",
            "3R\t250\t300\t3.7\n",
            "3R\t300\t400\t2.53\n",
            "3R\t400\t450\t4.28\n",
            "3R\t450\t500\t4.48\n",
            "3R\t500\t550\t2.14\n",
            "3R\t550\t600\t2.92\n",
            "3R\t600\t650\t5.45\n",
            "3R\t650\t750\t3.89\n",
            "3R\t750\t800\t2.73\n",
            "3R\t800\t900\t0\n",
            "3R\t900\t950\t0.58\n",
            "3R\t950\t1000\t1.36\n",
            "3R\t1000\t1050\t0.78\n",
            "3R\t1050\t1500\t0\n"
        ]

        assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
        unlink(outfile)


def test_bam_coverage_offset1():
    """
    Test -bs 1 --Offset 1
    """
    outfile = '/tmp/test_offset.bw'
    #for fname in [BAMFILE_A, CRAMFILE_A]:
    for fname in [BAMFILE_A]:
        args = "--Offset 1 --bam {} -p 1 -bs 1 -o {} -of bedgraph --verbose ".format(fname, outfile)
        print(args)
        args = args.split()
        bam_cov.main(args)
        _foo = open(outfile, 'r')
        resp = _foo.readlines()
        _foo.close()

        expected = [
            "3R\t0\t100\t0\n",
            "3R\t100\t101\t1\n",
            "3R\t101\t199\t0\n",
            "3R\t199\t200\t1\n",
            "chr_cigar\t0\t10\t0\n",
            "chr_cigar\t10\t11\t1\n",
            "chr_cigar\t11\t200\t0\n"
        ]

        assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
        unlink(outfile)


def test_bam_coverage_offset1_10():
    """
    Test -bs 1 --Offset 1 10
    """
    outfile = '/tmp/test_offset.bw'
    #for fname in [BAMFILE_A, CRAMFILE_A]:
    for fname in [BAMFILE_A]:
        args = "--Offset 1 10 -b {} -p 1 -bs 1 -of bedgraph -o {}".format(fname, outfile)
        args = args.split()
        bam_cov.main(args)
        _foo = open(outfile, 'r')
        resp = _foo.readlines()
        _foo.close()

        expected = [
            "3R\t0\t100\t0\n",
            "3R\t100\t110\t1\n",
            "3R\t110\t190\t0\n",
            "3R\t190\t200\t1\n",
            "chr_cigar\t0\t10\t0\n",
            "chr_cigar\t10\t20\t1\n",
            "chr_cigar\t20\t200\t0\n"
        ]

        assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
        unlink(outfile)


def test_bam_coverage_offset_minus1():
    """
    Test -bs 1 --Offset -1
    """
    outfile = '/tmp/test_offset.bw'
    #for fname in [BAMFILE_A, CRAMFILE_A]:
    for fname in [BAMFILE_A]:
        args = "--Offset -1 -b {} -p 1 -bs 1 -of bedgraph -o {}".format(fname, outfile)
        args = args.split()
        bam_cov.main(args)
        _foo = open(outfile, 'r')
        resp = _foo.readlines()
        _foo.close()

        expected = [
            "3R\t0\t149\t0\n",
            "3R\t149\t151\t1\n",
            "3R\t151\t200\t0\n",
            "chr_cigar\t0\t49\t0\n",
            "chr_cigar\t49\t50\t1\n",
            "chr_cigar\t50\t200\t0\n"
        ]

        assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
        unlink(outfile)


def test_bam_coverage_offset20_minus4():
    """
    Test -bs 1 --Offset 20 -4
    """
    outfile = '/tmp/test_offset.bw'
    #for fname in [BAMFILE_A, CRAMFILE_A]:
    for fname in [BAMFILE_A]:
        args = "--Offset 20 -4 -b {} -p 1 -bs 1 -of bedgraph -o {}".format(fname, outfile)
        args = args.split()
        bam_cov.main(args)
        _foo = open(outfile, 'r')
        resp = _foo.readlines()
        _foo.close()

        expected = [
            "3R\t0\t119\t0\n",
            "3R\t119\t147\t1\n",
            "3R\t147\t153\t0\n",
            "3R\t153\t181\t1\n",
            "3R\t181\t200\t0\n",
            "chr_cigar\t0\t29\t0\n",
            "chr_cigar\t29\t30\t1\n",
            "chr_cigar\t30\t40\t0\n",
            "chr_cigar\t40\t47\t1\n",
            "chr_cigar\t47\t200\t0\n"
        ]

        assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
        unlink(outfile)


def test_bam_compare_filter_blacklist():
    """
    Test --samFlagInclude --samFlagExclude --minMappingQuality --ignoreDuplicates and --blackListFileName
    """
    outfile = '/tmp/test_file_filter.bg'
    #for A, B in [(BAMFILE_FILTER1, BAMFILE_FILTER2), (CRAMFILE_FILTER1, CRAMFILE_FILTER2)]:
    for A, B in [(BAMFILE_FILTER1, BAMFILE_FILTER2)]:
        args = "-b1 {} -b2 {} -p 1 -o {} -of bedgraph --samFlagInclude 512 " \
               "--samFlagExclude 256 --minMappingQuality 5 --verbose " \
               "--blackListFileName {}".format(A, B, outfile, BEDFILE_FILTER)
        print(args)
        args = args.split()
        bam_comp.main(args)

        _foo = open(outfile, 'r')
        resp = _foo.readlines()
        _foo.close()
        expected = [
            "3R\t0\t100\t0\n",
            "3R\t100\t150\t-0.13\n",
            "3R\t150\t200\t-0.1\n",
            "3R\t200\t250\t0.03\n",
            "3R\t250\t300\t0.21\n",
            "3R\t300\t350\t0.09\n",
            "3R\t350\t400\t-0.1\n",
            "3R\t400\t450\t-0.01\n",
            "3R\t450\t500\t0.05\n",
            "3R\t500\t550\t0.11\n",
            "3R\t550\t600\t0.08\n",
            "3R\t600\t700\t-0.01\n",
            "3R\t700\t750\t-0.13\n",
            "3R\t750\t800\t-0.18\n",
            "3R\t800\t900\t0\n",
            "3R\t900\t950\t0.32\n",
            "3R\t950\t1000\t0.16\n",
            "3R\t1000\t1050\t-0.01\n",
            "3R\t1050\t1500\t0\n"
        ]
        assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
        unlink(outfile)
