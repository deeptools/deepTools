import deeptools.bamCoverage as bam_cov
import deeptools.bamCompare as bam_comp
from os import unlink

ROOT = "./test/test_data/"
BAMFILE_A  = ROOT + "testA.bam"
BAMFILE_B  = ROOT + "testB.bam"


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

    resp =open(outfile, 'r').readlines()
    assert resp == ['3R\t0\t50\t0.00\n', '3R\t50\t150\t1.00\n', '3R\t150\t200\t2.0\n']
    unlink(outfile)

def test_bam_compare_arguments():
    """
    Test minimal command line args for bamCoverage
    """
    outfile = '/tmp/test_file.bg'
    args = "--bamfile1 {} --bamfile2 {} " \
           "-o {} --outFileFormat bedgraph".format(BAMFILE_A, BAMFILE_B, outfile).split()
    bam_comp.main(args)

    resp =open(outfile, 'r').readlines()
    print resp
    expected = ['3R\t0\t50\t0.00\n', '3R\t50\t100\t-0.81\n', '3R\t100\t150\t0.19\n', '3R\t150\t200\t-0.3\n']
    assert resp == expected
    unlink(outfile)


def test_ignore_for_normalization():
    """
    test does not work because a weird behaviour of pysam.idxstats under nosetest
    """
    pass
    #outfile = '/tmp/test_file.bg'
    #args = "--bam {} -o {} --outFileFormat bedgraph " \
    #       "--ignoreForNormalization chr_cigar --normalizeUsingRPKM".format(BAMFILE_A, outfile).split()
    #bam_cov.main(args)

    #resp =open(outfile, 'r').readlines()
    #assert resp == ['3R\t0\t100\t0.00\n', '3R\t100\t200\t1.0\n', 'chr_cigar\t0\t50\t2.00\n']
