import os
import filecmp
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile

import deeptools.plotEnrichment
import deeptools.utilities


__author__ = 'Bjoern'


TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_plotEnrichment/"

# Headroom for anti-aliasing/text-hinting drift between matplotlib patch
# versions (e.g. 3.10.x vs 3.11.x render identical plots with slightly
# different sub-pixel edges); see matplotlib_defaults.py.
tolerance = 20


def run_plotEnrichment(extra_args=None):
    """Run plotEnrichment and return generated plot and count files."""

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
        f"--BED {TEST_DATA}test.gtf "
        f"--plotFile {plotfile.name} "
        f"--plotFileFormat png "
        f"--outRawCounts {txtfile.name}"
    ).split()

    if extra_args:
        args.extend(extra_args)

    deeptools.plotEnrichment.main(args)

    return plotfile.name, txtfile.name


def cleanup(*files):
    for f in files:
        if os.path.exists(f):
            os.remove(f)


def test_plotEnrichment_default():
    """Test default plotEnrichment output and raw counts."""

    plotfile, txtfile = run_plotEnrichment()

    try:
        assert filecmp.cmp(
            os.path.join(ROOT, "outRawCounts_default.tabular"),
            txtfile
        ) is True

        res = compare_images(
            os.path.join(ROOT, "plotEnrichment_defaults.png"),
            plotfile,
            tolerance
        )

        assert res is None, res

    finally:
        cleanup(plotfile, txtfile)


def test_plotEnrichment_ggplot():
    """Test plotEnrichment output with --ggplot."""

    plotfile, txtfile = run_plotEnrichment(
        ["--ggplot"]
    )

    try:
        res = compare_images(
            os.path.join(ROOT, "plotEnrichment_ggplot.png"),
            plotfile,
            tolerance
        )

        assert res is None, res

    finally:
        cleanup(plotfile, txtfile)