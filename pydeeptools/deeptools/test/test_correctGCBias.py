# from unittest import TestCase

import deeptools.correctGCBias
import deeptools.utilities
import os.path
from os import unlink
import pysam

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"


def test_correctGCBias():
    """
    Test minimal command line args for correctGCBias
    """
    GCbiasFreq = ROOT + 'computeGCBias_result1.tabular'
    BAM = ROOT + 'paired_chr2L.bam'
    GENOME = ROOT + 'sequence.2bit'
    outfile = '/tmp/test_correctbias.bam'
    args = "--GCbiasFrequenciesFile {} --bamfile {} --genome {} --effectiveGenomeSize 10050 --correctedFile {}".format(GCbiasFreq, BAM, GENOME, outfile).split()
    deeptools.correctGCBias.main(args)

    alignment_count = pysam.AlignmentFile(outfile, "rb").count()
    expected_count = 11630
    assert abs(alignment_count - expected_count) < 50
    unlink(outfile)