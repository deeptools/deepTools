import os
import filecmp
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import deeptools.plotCoverage
import deeptools.utilities

__author__ = 'Bjoern'

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_plotCoverage/"

tolerance = 13  # default matplotlib pixed difference tolerance


def test_plotCoverage_default():
    plotfile = NamedTemporaryFile(suffix='.png', prefix='deeptools_testfile_', delete=False)
    txtfile = NamedTemporaryFile(suffix='.tab', prefix='deeptools_testfile_', delete=False)

    args = "--bamfiles {0}test1.bam {0}test2.bam --plotFile {1}" \
           " --plotFileFormat png --outRawCounts {2}".format(TEST_DATA, plotfile.name, txtfile.name).split()
    deeptools.plotCoverage.main(args)
    assert filecmp.cmp(os.path.join(ROOT, 'outRawCounts_default.tabular'), txtfile.name) is True

    res = compare_images(ROOT + 'plotCoverage_default.png', plotfile.name, tolerance)
    assert res is None, res
    os.remove(txtfile.name)
    os.remove(plotfile.name)
