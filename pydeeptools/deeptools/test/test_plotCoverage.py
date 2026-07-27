import os
import filecmp
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile

import deeptools.plotCoverage
import deeptools.utilities


__author__ = 'Bjoern'


TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_plotCoverage/"

tolerance = 13


def run_plotCoverage(extra_args=None):
    """Run plotCoverage and return generated plot and raw count files."""

    plotfile = NamedTemporaryFile(
        suffix='.png',
        prefix='deeptools_testfile_',
        delete=False
    )

    txtfile = NamedTemporaryFile(
        suffix='.tab',
        prefix='deeptools_testfile_',
        delete=False
    )

    args = (
        f"--bamfiles {TEST_DATA}test1.bam "
        f"{TEST_DATA}test2.bam "
        f"--plotFile {plotfile.name} "
        f"--plotFileFormat png "
        f"--outRawCounts {txtfile.name}"
    ).split()

    if extra_args:
        args.extend(extra_args)

    deeptools.plotCoverage.main(args)

    return plotfile.name, txtfile.name


def cleanup(*files):
    for f in files:
        if os.path.exists(f):
            os.remove(f)


def test_plotCoverage_default():
    """Test default plotCoverage output and raw counts."""

    plotfile, txtfile = run_plotCoverage()

    try:
        assert filecmp.cmp(
            os.path.join(ROOT, "outRawCounts_default.tabular"),
            txtfile
        ) is True

        res = compare_images(
            os.path.join(ROOT, "plotCoverage_default.png"),
            plotfile,
            tolerance
        )

        assert res is None, res

    finally:
        cleanup(plotfile, txtfile)


def test_plotCoverage_ggplot():
    """Test plotCoverage output with --ggplot."""

    plotfile, txtfile = run_plotCoverage(
        ["--ggplot"]
    )

    try:
        res = compare_images(
            os.path.join(ROOT, "plotCoverage_ggplot.png"),
            plotfile,
            tolerance
        )

        assert res is None, res

    finally:
        cleanup(plotfile, txtfile)