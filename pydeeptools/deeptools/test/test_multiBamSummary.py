import deeptools.multiBamSummary2 as mbs
import numpy as np
import numpy.testing as nt

import os.path
from os import unlink
import tempfile

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/test_mbs/"
BAM = ROOT + "test1"
GTF = ROOT + "test.gtf"
BED = ROOT + "test.bed"
BAMA = ROOT + "testA"
BAMB = ROOT + "testB"



def test_multiBamSummary_bedmode_gtf():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"BED-file --BED {GTF} -b {fname} {fname} -o {outfile}".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix, np.array([[144.0, 144.0],[143.0, 143.0]]))

def test_multiBamSummary_bedmode_bed():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"BED-file --BED {BED} -b {fname} {fname} -o {outfile}".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix,
            np.array(
                [[1.0, 1.0], [144.0, 144.0], [144.0, 144.0], [6.0, 6.0], [143.0, 143.0], [22.0, 22.0], [25.0, 25.0], [1.0, 1.0], [0.0, 0.0]]
            )
        )

def test_multiBamSummary_bedmode_multibed():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"BED-file --BED {BED} {GTF} {BED} -b {fname} {fname} -o {outfile}".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix,
            np.array(
                [[144.0, 144.0], [1.0, 1.0], [1.0, 1.0], [144.0, 144.0], [144.0, 144.0], [144.0, 144.0], [144.0, 144.0], [143.0, 143.0], [6.0, 6.0], [6.0, 6.0], [143.0, 143.0], [143.0, 143.0], [22.0, 22.0], [22.0, 22.0], [25.0, 25.0], [25.0, 25.0], [1.0, 1.0], [1.0, 1.0], [0.0, 0.0], [0.0, 0.0]]
            )
        )


def test_multiBamSummary_metagene():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f'BED-file --BED {GTF} -b {fname} {fname} -o {outfile} --metagene'.split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix, np.array([[25.0, 25.0],
                                             [31.0, 31.0]]))
        unlink(outfile)


def test_multiBamSummary_scalingFactors():
    _, outfile = tempfile.mkstemp(suffix=".txt")
    _, outfile2 = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        bama = BAMA + fname
        bamb = BAMB + fname
        args = f'bins --binSize 50 -b {bama} {bamb} --scalingFactors {outfile} -o {outfile2} --verbose'.split()
        mbs.main(args)
        resp = open(outfile).read().strip().split('\n')
        assert resp == ["Sample\tscalingFactor", f"testA{fname}\t1.1892071", f"testB{fname}\t0.8408964"]
        nt.assert_equal(resp, ["Sample\tscalingFactor", f"testA{fname}\t1.1892071", f"testB{fname}\t0.8408964"])
