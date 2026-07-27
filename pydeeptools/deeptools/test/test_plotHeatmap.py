import os
import tempfile

import matplotlib.pyplot as plt
from matplotlib.testing.compare import compare_images

import deeptools.plotHeatmap


TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data")
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_plotHeatmap")

tolerance = 13


def run_plotHeatmap(args):
    """Run plotHeatmap and return generated output file."""

    fd, plotfile = tempfile.mkstemp(suffix=".png")
    os.close(fd)

    args = args + [
        "--outFileName",
        plotfile
    ]

    plt.ioff()
    deeptools.plotHeatmap.main(args)
    
    return plotfile


def cleanup(plotfile):
    plt.close("all")
    if os.path.exists(plotfile):
        os.remove(plotfile)


def test_plotHeatmap_default():
    """Image comparison test for default output."""

    args = (
        f"--matrixFile {TEST_DATA}/computeMatrix_result1.gz"
    ).split()

    plotfile = run_plotHeatmap(args)

    try:
        res = compare_images(
            os.path.join(ROOT, "plotHeatmap_default.png"),
            plotfile,
            tolerance
        )

        assert res is None, (
            f"{plotfile} doesn't match "
            f"{ROOT}/plotHeatmap_default.png"
        )

    finally:
        cleanup(plotfile)


def test_plotHeatmap_ggplot():
    """Check that --ggplot produces an output image."""

    args = (
        f"--matrixFile {TEST_DATA}/computeMatrix_result1.gz "
        "--ggplot"
    ).split()

    plotfile = run_plotHeatmap(args)

    try:
        res = compare_images(
            os.path.join(ROOT, "plotHeatmap_ggplot.png"),
            plotfile,
            tolerance
        )

        assert res is None, (
            f"{plotfile} doesn't match "
            f"{ROOT}/plotHeatmap_ggplot.png"
        )


    finally:
        cleanup(plotfile)