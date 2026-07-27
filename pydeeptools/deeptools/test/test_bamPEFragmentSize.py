import os.path
import tempfile

from matplotlib.testing.compare import compare_images

import deeptools.bamPEFragmentSize


ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data"

tolerance = 10


def run_bamPEFragmentSize_histogram(extra_args=None):
    """
    Run bamPEFragmentSize histogram and return generated png file.
    """

    fd, outfile = tempfile.mkstemp(suffix=".png")
    os.close(fd)

    args = (
        f"--bamfiles {ROOT}/bowtie2_test1.bam "
        f"--samplesLabel bowtie2_test1.bam "
        f"--plotFileFormat png "
        f"--plotTitle Test-Plot "
        f"--histogram {outfile}"
    ).split()

    if extra_args:
        args.extend(extra_args)

    deeptools.bamPEFragmentSize.main(args)

    return outfile


def cleanup(*files):
    """
    Remove temporary files.
    """
    for f in files:
        if f and os.path.exists(f):
            os.remove(f)


def test_bamPEFragmentSize_histogram():
    """
    Test histogram plot for bamPEFragmentSize.
    """

    outfile = run_bamPEFragmentSize_histogram()

    try:
        res = compare_images(
            ROOT + "/bamPEFragmentSize_histogram_result1.png",
            outfile,
            tolerance
        )

        assert res is None, res

    finally:
        cleanup(outfile)


def test_bamPEFragmentSize_histogram_ggplot():
    """
    Test histogram plot for bamPEFragmentSize using --ggplot.
    """

    outfile = run_bamPEFragmentSize_histogram(
        ["--ggplot"]
    )

    try:
        res = compare_images(
            ROOT + "/bamPEFragmentSize_histogram_ggplot.png",
            outfile,
            tolerance
        )

        assert res is None, res

    finally:
        print(outfile)


def test_bamPEFragmentSize_fr_sizes():
    """
    Test fragment length information for bamPEFragmentSize.
    """

    fd, out_lengths = tempfile.mkstemp(suffix=".txt")
    os.close(fd)

    fd, out_metrics = tempfile.mkstemp(suffix=".txt")
    os.close(fd)

    args = (
        f"--bamfiles {ROOT}/bowtie2_test1.bam "
        f"--outRawFragmentLengths {out_lengths} "
        f"--table {out_metrics}"
    ).split()

    deeptools.bamPEFragmentSize.main(args)

    try:

        with open(out_lengths, "r") as l:
            l_resp = l.readlines()

        l_expected = [
            "241\t1",
            "242\t1",
            "251\t1"
        ]

        matches = [
            expected
            for expected in l_expected
            if any(expected in resp for resp in l_resp)
        ]

        assert matches == l_expected


        with open(out_metrics, "r") as m:
            m_resp = m.readlines()[1]


        m_expected = (
            "3\t241.0\t241.5\t244.66666666666666\t242.0\t246.5\t251.0\t"
            "4.496912521077347\t1.0\t241.2\t241.4\t241.6\t241.8\t243.8\t245.6\t"
            "247.4\t249.2\t250.82\t3\t251.0\t251.0\t251.0\t251.0\t251.0\t251.0\t"
            "0.0\t0.0\t251.0\t251.0\t251.0\t251.0\t251.0\t251.0\t251.0\t251.0\t251.0\n"
        )

        assert m_expected in m_resp

    finally:
        cleanup(
            out_lengths,
            out_metrics
        )