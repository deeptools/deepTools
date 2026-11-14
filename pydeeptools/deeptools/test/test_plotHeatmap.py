import os
import tempfile

import matplotlib.pyplot as plt
from matplotlib.testing.compare import compare_images

import deeptools.plotHeatmap

TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data")
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_plotHeatmap")

# Headroom for anti-aliasing/text-hinting drift between matplotlib patch
# versions (e.g. 3.10.x vs 3.11.x render identical plots with slightly
# different sub-pixel edges); see matplotlib_defaults.py.
tolerance = 25


def run_plotHeatmap(args):
    """Run plotHeatmap and return generated output file."""

    _fd, plotfile = tempfile.mkstemp(suffix=".png")
    # os.close(fd)

    args = args + [
        "--outFileName",
        plotfile
    ]

    plt.ioff()
    deeptools.plotHeatmap.main(args)

    return plotfile


# def run_plotHeatmap(args, Typ):

#     plotfile = os.path.join(
#         ROOT,
#         "debug"+ Typ + ".png"
#     )

#     if os.path.exists(plotfile):
#         os.remove(plotfile)

#     args = args + [
#         "--outFileName",
#         plotfile
#     ]

#     plt.ioff()
#     deeptools.plotHeatmap.main(args)

#     print("Generated:", plotfile)

#     return plotfile


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

def _run_silhouette(method):
    """Run plotHeatmap with clustering + --silhouette and return the sorted-regions BED header."""
    _fd, bedfile = tempfile.mkstemp(suffix=".bed")
    args = (
        f"--matrixFile {TEST_DATA}/computeMatrix_result1.gz "
        f"{method} 2 --silhouette --outFileSortedRegions {bedfile}"
    ).split()

    plotfile = run_plotHeatmap(args)
    try:
        with open(bedfile) as fh:
            return fh.readline().rstrip("\n").split("\t")
    finally:
        cleanup(plotfile)
        os.remove(bedfile)


def test_plotHeatmap_silhouette_kmeans():
    assert _run_silhouette("--kmeans")[-1] == "silhouette"


def test_plotHeatmap_silhouette_hclust():
    assert _run_silhouette("--hclust")[-1] == "silhouette"
