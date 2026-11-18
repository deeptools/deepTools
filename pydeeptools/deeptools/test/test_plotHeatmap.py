import os
import tempfile

import matplotlib.pyplot as plt
import pyBigWig
import pytest
from matplotlib.testing.compare import compare_images

import deeptools.computeMatrix2
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


def test_plotHeatmap_silhouette_two_regions_and_singleton(tmp_path):
    track = tmp_path / "signal.bw"
    with pyBigWig.open(str(track), "w") as bw:
        bw.addHeader([("chr1", 100)])
        bw.addEntries(["chr1"] * 3, [0, 30, 60],
                      ends=[20, 50, 80], values=[1.0, 2.0, 10.0])
    regions = tmp_path / "regions.bed"
    regions.write_text("chr1\t0\t20\tlow\t0\t+\nchr1\t30\t50\tmid\t0\t+\nchr1\t60\t80\thigh\t0\t+\n")
    matrix = tmp_path / "matrix.gz"
    deeptools.computeMatrix2.main([
        "reference-point", "-S", str(track), "-R", str(regions), "-b", "0",
        "-a", "20", "-bs", "10", "-o", str(matrix),
    ])
    sorted_regions = tmp_path / "sorted.bed"
    deeptools.plotHeatmap.main([
        "-m", str(matrix), "-out", str(tmp_path / "heatmap.png"), "--hclust", "2",
        "--silhouette", "--outFileSortedRegions", str(sorted_regions),
    ])

    rows = [line.split("\t") for line in sorted_regions.read_text().splitlines()[1:]]
    groups = {row[3]: row[-2] for row in rows}
    assert groups["low"] == groups["mid"] != groups["high"]
    scores = {row[3]: float(row[-1]) for row in rows}
    assert scores == pytest.approx({"low": 8 / 9, "mid": 7 / 8, "high": 0.0})
