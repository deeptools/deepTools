import os
from tempfile import NamedTemporaryFile

from matplotlib.testing.compare import compare_images

import deeptools.plotProfile


TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_plotProfile/"

tolerance = 13


def run_plotProfile(args):
    """Run plotProfile and return generated plot file."""

    plotfile = NamedTemporaryFile(
        suffix='.png',
        prefix='deeptools_testfile_',
        delete=False
    )

    args += [
        "--outFileName",
        plotfile.name,
        "--plotFileFormat",
        "png"
    ]

    deeptools.plotProfile.main(args)

    return plotfile.name


def cleanup(plotfile):
    if os.path.exists(plotfile):
        os.remove(plotfile)


def test_plotProfile_default():
    """Image comparison test for default plotProfile output."""

    args = (
        f"--matrixFile {TEST_DATA}computeMatrix_result1.gz"
    ).split()

    plotfile = run_plotProfile(args)

    try:
        res = compare_images(
            ROOT + "plotProfile_default.png",
            plotfile,
            tolerance
        )

        assert res is None, res

    finally:
        print(plotfile)


def test_plotProfile_ggplot():
    """Image comparison test for --ggplot output."""

    args = (
        f"--matrixFile {TEST_DATA}computeMatrix_result1.gz "
        "--ggplot"
    ).split()

    plotfile = run_plotProfile(args)

    try:
        res = compare_images(
            ROOT + "plotProfile_ggplot.png",
            plotfile,
            tolerance
        )

        assert res is None, res

    finally:
        print(plotfile)