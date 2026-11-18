import os
import tempfile

import matplotlib.pyplot as plt
import pytest
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

def test_plotHeatmap_clustering_preserves_missing_signal(tmp_path):
    import gzip

    import numpy as np

    matrix = os.path.join(TEST_DATA, "computeMatrixOperations.mat.gz")
    saved = tmp_path / "clustered.gz"
    deeptools.plotHeatmap.main([
        "-m", matrix, "-out", str(tmp_path / "heatmap.png"), "--hclust", "2",
        "--outFileNameMatrix", str(saved),
    ])

    def read_rows(path):
        with gzip.open(path, "rt") as handle:
            next(handle)
            rows = (line.rstrip("\n").split("\t") for line in handle)
            return {tuple(row[:4]): np.array(row[6:], dtype=float) for row in rows}

    original = read_rows(matrix)
    clustered = read_rows(saved)
    assert len(original) == 196
    assert sum(np.isnan(row).sum() for row in original.values()) == 4822
    assert clustered.keys() == original.keys()
    for region, values in original.items():
        np.testing.assert_array_equal(clustered[region], values)


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


def test_plotHeatmap_colors():
    """--colors sets the profile line colors (recycled) and produces a plot."""
    args = (
        f"--matrixFile {TEST_DATA}/computeMatrix_result1.gz "
        "--colors red #00ff00"
    ).split()
    plotfile = run_plotHeatmap(args)
    try:
        assert os.path.getsize(plotfile) > 0
    finally:
        cleanup(plotfile)


def test_plotHeatmap_colors_invalid():
    """An invalid --colors value is rejected by argparse."""
    _fd, plotfile = tempfile.mkstemp(suffix=".png")
    os.close(_fd)
    args = (
        f"--matrixFile {TEST_DATA}/computeMatrix_result1.gz "
        f"--colors notacolor --outFileName {plotfile}"
    ).split()
    try:
        with pytest.raises(SystemExit):
            deeptools.plotHeatmap.main(args)
    finally:
        cleanup(plotfile)
