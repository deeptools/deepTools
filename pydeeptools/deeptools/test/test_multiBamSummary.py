import deeptools.multiBamSummary2 as mbs
import numpy as np
import numpy.testing as nt

import os.path
from os import unlink
import tempfile

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
BAM = ROOT + "test1.bam"
CRAM = ROOT + "test1.cram"
GTF = ROOT + "test.gtf"
BAMA = ROOT + "testA.bam"
BAMB = ROOT + "testB.bam"


def test_multiBamSummary_gtf():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    #for fname in [BAM, CRAM]:
    for fname in [BAM]:
        args = 'BED-file --BED {0} -b {1} {1} -o {2}'.format(GTF, fname, outfile).split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix, np.array([[144.0, 144.0],
                                             [143.0, 143.0]]))
        unlink(outfile)


def test_multiBamSummary_metagene():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    #for fname in [BAM, CRAM]:
    for fname in [BAM]:
        args = 'BED-file --BED {0} -b {1} {1} -o {2} --metagene'.format(GTF, fname, outfile).split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix, np.array([[24.0, 24.0],
                                             [31.0, 31.0]]))
        unlink(outfile)


def test_multiBamSummary_scalingFactors():
    _, outfile = tempfile.mkstemp(suffix=".txt")
    _, outfile2 = tempfile.mkstemp(suffix=".npz")
    args = 'bins --binSize 50 -b {} {} --scalingFactors {} -o {} --verbose'.format(BAMA, BAMB, outfile, outfile2).split()
    mbs.main(args)
    resp = open(outfile).read().strip().split('\n')
    assert resp == ["Sample\tscalingFactor", "testA.bam\t1.1892071", "testB.bam\t0.8408964"]
    nt.assert_equal(resp, ["Sample\tscalingFactor", "testA.bam\t1.1892071", "testB.bam\t0.8408964"])
    unlink(outfile)
    unlink(outfile2)
