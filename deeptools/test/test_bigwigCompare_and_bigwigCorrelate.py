import deeptools.bigwigCompare as bwComp
import os.path
from os import unlink

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
BIGWIG_A = ROOT + "testA_skipNAs.bw"
BIGWIG_B = ROOT + "testB_skipNAs.bw"


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



def test_bigwigCompare():
    outfile = '/tmp/result.bg'
    args = "-b1 {} -b2 {} -o {} --ratio add --outFileFormat bedgraph".format(BIGWIG_A, BIGWIG_B, outfile).split()
    print " ".join(args)
    bwComp.main(args)
    resp = open(outfile, 'r').readlines()
    expected = ['3R\t0\t50\t0.00\n', '3R\t50\t100\t1.00\n', '3R\t100\t150\t2.00\n', '3R\t150\t200\t3.0\n']
    assert resp == expected, "{} != {}".format(resp, expected)
    unlink(outfile)

def test_bigwigCompare_skipnas():
    outfile = '/tmp/result.bg'
    args = "-b1 {} -b2 {} -o {} --ratio add --skipNAs " \
           "--outFileFormat bedgraph".format(BIGWIG_A, BIGWIG_B, outfile).split()
    print " ".join(args)
    bwComp.main(args)
    resp = open(outfile, 'r').readlines()
    expected = ['3R\t100\t150\t2.00\n', '3R\t150\t200\t3.0\n']
    assert resp == expected, "{} != {}".format(resp, expected)
    unlink(outfile)

#class TestBigwigTools(object):

#    def setUp(self):
#        # create bigwig files based on bam files
#        outfile = '/tmp/test_file.bw'
#        args = "--bam {} -o {} --outFileFormat bedgraph".format(BAMFILE_B, outfile).split()
#        bam_cov.main(args)
