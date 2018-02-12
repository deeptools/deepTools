import deeptools.multiBamSummary as mbs
import numpy as np
import numpy.testing as nt

import os.path
from os import unlink

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
BAM = ROOT + "test1.bam"
CRAM = ROOT + "test1.cram"
GTF = ROOT + "test.gtf"


def test_multiBamSummary_gtf():
    outfile = '/tmp/_test.npz'
    for fname in [BAM, CRAM]:
        args = 'BED-file --BED {0} -b {1} {1} -o {2}'.format(GTF, fname, outfile).split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']
        labels = resp['labels']
        if fname == BAM:
            nt.assert_equal(labels, ['test1.bam', 'test1.bam'])
        else:
            nt.assert_equal(labels, ['test1.cram', 'test1.cram'])
        nt.assert_allclose(matrix, np.array([[144.0, 144.0],
                                             [143.0, 143.0]]))
        unlink(outfile)


def test_multiBamSummary_metagene():
    outfile = '/tmp/_test.npz'
    for fname in [BAM, CRAM]:
        args = 'BED-file --BED {0} -b {1} {1} -o {2} --metagene'.format(GTF, fname, outfile).split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']
        labels = resp['labels']
        if fname == BAM:
            nt.assert_equal(labels, ['test1.bam', 'test1.bam'])
        else:
            nt.assert_equal(labels, ['test1.cram', 'test1.cram'])
        nt.assert_allclose(matrix, np.array([[25.0, 25.0],
                                             [31.0, 31.0]]))
        unlink(outfile)
