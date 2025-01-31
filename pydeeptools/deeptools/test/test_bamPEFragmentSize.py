import deeptools.bamPEFragmentSize
from matplotlib.testing.compare import compare_images
import os.path
import filecmp
from os import unlink
from tempfile import NamedTemporaryFile

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data"


def test_bamPEFragmentSize_histogram():
    """
    Test histogram plot for bamPEFragmentSize
    """
    outfile = '/tmp/test_histogram.png'
    args = "--bamfiles {}/bowtie2_test1.bam --samplesLabel bowtie2_test1.bam --plotFileFormat png --plotTitle Test-Plot --histogram {}".format(ROOT, outfile).split()
    deeptools.bamPEFragmentSize.main(args)

    res = compare_images(ROOT + '/bamPEFragmentSize_histogram_result1.png', outfile, 10)
    assert res is None, res
    unlink(outfile)

def test_bamPEFragmentSize_fr_sizes():
    """
    Test fragment length information for bamPEFragmentSize
    """
    out_lengths = '/tmp/test_raw_frag_lengths.txt'
    out_metrics = '/tmp/test_metrics_table.txt'
    args = "--bamfiles {}/bowtie2_test1.bam --outRawFragmentLengths {} --table {}".format(ROOT, out_lengths, out_metrics).split()
    deeptools.bamPEFragmentSize.main(args)

    l = open(out_lengths, 'r')
    l_resp = l.readlines()
    l.close()
    l_expected = ['241\t1', '242\t1', '251\t1']
    matches = [expected for expected in l_expected if any(expected in resp for resp in l_resp)]
    assert matches == l_expected

    m = open(out_metrics, 'r')
    m_resp = m.readlines()[1]
    m.close()
    m_expected = '3\t241.0\t241.5\t244.66666666666666\t242.0\t246.5\t251.0\t4.496912521077347\t1.0\t241.2\t241.4\t241.6\t241.8\t243.8\t245.6\t247.4\t249.2\t250.82\t3\t251.0\t251.0\t251.0\t251.0\t251.0\t251.0\t0.0\t0.0\t251.0\t251.0\t251.0\t251.0\t251.0\t251.0\t251.0\t251.0\t251.0\n'
    assert m_expected in f"{m_resp}"
    unlink(out_lengths)
    unlink(out_metrics)