import deeptools.bigwigAverage as bwAve

import os.path
from os import unlink

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
BIGWIG_A = ROOT + "testA_skipNAs.bw"
BIGWIG_B = ROOT + "testB_skipNAs.bw"
BIGWIG_C = ROOT + "test1.bw.bw"


"""
The distribution of reads for the bam file is:

              0                              100                           200
              |------------------------------------------------------------|
testA.bam  3R                                ==============>
                                                            <==============


testB.bam  3R                 <==============               ==============>
                                             ==============>
                                                            ==============>

The resulting bigwig files are as follows:

testA_skipNas:
    3R      100     200     1
    chr_cigar       0       50      2

testB_skipNas:
    3R      50      150     1
    3R      150     200     2
"""


def test_bigwigAverage():
    outfile = '/tmp/result.bg'
    args = "--bigwigs {} {} -o {} --outFileFormat bedgraph".format(BIGWIG_A, BIGWIG_B, outfile).split()
    bwAve.main(args)
    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t0\t50\t0\n', '3R\t50\t100\t0.5\n', '3R\t100\t150\t1\n', '3R\t150\t200\t1.5\n']
    assert resp == expected, "{} != {}".format(resp, expected)
    unlink(outfile)


def test_bigwigAverage_skipnas():
    outfile = '/tmp/result.bg'
    args = "--bigwigs {} {} -o {} --skipNAs " \
           "--outFileFormat bedgraph".format(BIGWIG_A, BIGWIG_B, outfile).split()
    bwAve.main(args)
    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t100\t150\t1\n', '3R\t150\t200\t1.5\n']
    assert resp == expected, "{} != {}".format(resp, expected)
    unlink(outfile)


def test_bigwigAverageWithScale():
    outfile = '/tmp/result.bg'
    args = "--bigwigs {} {} -o {} --outFileFormat bedgraph --scaleFactors 1:0.5".format(BIGWIG_A, BIGWIG_B, outfile).split()
    bwAve.main(args)
    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['3R\t0\t50\t0\n', '3R\t50\t100\t0.25\n', '3R\t100\t150\t0.75\n', '3R\t150\t200\t1\n']
    assert resp == expected, "{} != {}".format(resp, expected)
    unlink(outfile)
