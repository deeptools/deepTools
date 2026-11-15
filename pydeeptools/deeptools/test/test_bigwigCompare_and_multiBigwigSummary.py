import os.path
import tempfile
from os import unlink

import numpy as np
import numpy.testing as nt

import deeptools.bigwigCompare as bwComp
import deeptools.multiBigwigSummary as bwCorr

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


def test_bigwigCompare():
    _, outfile = tempfile.mkstemp(suffix=".bg")
    args = f"-b1 {BIGWIG_A} -b2 {BIGWIG_B} -o {outfile} --operation add --outFileFormat bedgraph".split()
    bwComp.main(args)
    with open(outfile, 'r') as _foo:
        resp = _foo.readlines()
    expected = ['3R\t0\t50\t0\n', '3R\t50\t100\t1\n', '3R\t100\t150\t2\n', '3R\t150\t200\t3\n']
    assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
    unlink(outfile)


def test_bigwigCompare_skipnas():
    _, outfile = tempfile.mkstemp(suffix=".bg")
    args = f"-b1 {BIGWIG_A} -b2 {BIGWIG_B} -o {outfile} --operation add --skipNAs " \
           "--outFileFormat bedgraph".split()
    bwComp.main(args)
    with open(outfile, 'r') as _foo:
        resp = _foo.readlines()
    expected = ['3R\t100\t150\t2\n', '3R\t150\t200\t3\n']
    assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
    unlink(outfile)


def test_bigwigCompare_skipZeroOverZero():
    _, outfile = tempfile.mkstemp(suffix=".bg")
    args = f"-b1 {BIGWIG_A} -b2 {BIGWIG_A} -o {outfile} --skipZeroOverZero --pseudocount 1 3 --outFileFormat bedgraph".split()
    bwComp.main(args)
    with open(outfile, 'r') as _foo:
        resp = _foo.readlines()
    expected = ['3R\t100\t200\t-1\n']
    assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
    unlink(outfile)


def test_multiBigwigSummary():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    args = f"bins -b {BIGWIG_A} {BIGWIG_B} --binSize 50 -o {outfile}".split()
    bwCorr.main(args)
    resp = np.load(outfile)
    matrix = resp['matrix']
    labels = resp['labels']
    nt.assert_equal(matrix, np.array([[np.nan, np.nan],
                                      [np.nan, 1.],
                                      [1., 1.],
                                      [1., 2.]]))
    nt.assert_equal(labels, ['testA_skipNAs.bw', 'testB_skipNAs.bw'])
    unlink(outfile)


def test_multiBigwigSummary_outrawcounts():
    """
    Test multiBigwigSummary raw counts output
    """
    _, nullfile = tempfile.mkstemp(suffix=".npz")
    _, outfile = tempfile.mkstemp(suffix=".txt")
    args = f"bins -b {BIGWIG_A} {BIGWIG_B} --binSize 50 -o {nullfile} --outRawCounts {outfile} ".split()
    bwCorr.main(args)
    with open(outfile, 'r') as _foo:
        resp = _foo.read()
    expected = """#'chr'	'start'	'end'	'testA_skipNAs.bw'	'testB_skipNAs.bw'
3R	0	50	nan	nan
3R	50	100	nan	1.0
3R	100	150	1.0	1.0
3R	150	200	1.0	2.0
"""
    assert f"{resp}" == f"{expected}", f"{resp} != {expected}"
    unlink(outfile)
    unlink(nullfile)


def test_multiBigwigSummary_gtf():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    args = f"BED-file -b {BIGWIG_C} {BIGWIG_C} --BED {ROOT}/test.gtf -o {outfile}".split()
    bwCorr.main(args)
    resp = np.load(outfile)
    matrix = resp['matrix']
    labels = resp['labels']
    nt.assert_equal(labels, ['test1.bw.bw', 'test1.bw.bw'])
    nt.assert_allclose(matrix, np.array([[27.475, 27.475],
                                         [27.31248719, 27.31248719]]))
    unlink(outfile)


def test_multiBigwigSummary_metagene():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    args = f"BED-file --metagene -b {BIGWIG_C} {BIGWIG_C} --BED {ROOT}/test.gtf -o {outfile}".split()
    bwCorr.main(args)
    resp = np.load(outfile)
    matrix = resp['matrix']
    labels = resp['labels']
    nt.assert_equal(labels, ['test1.bw.bw', 'test1.bw.bw'])
    nt.assert_allclose(matrix, np.array([[20.28956028, 20.28956028],
                                         [22.1923501, 22.1923501]]))
    unlink(outfile)
