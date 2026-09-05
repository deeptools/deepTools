import os
from tempfile import NamedTemporaryFile

from matplotlib.testing.compare import compare_images

import deeptools.plotFingerprint


TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_plotFingerprint/"

tolerance = 13


def run_plotFingerprint(args):
    """Run plotFingerprint and return generated plot file."""

    plotfile = NamedTemporaryFile(
        suffix=".png",
        prefix="deeptools_testfile_",
        delete=False
    )

    args.extend([
        "-o",
        plotfile.name,
        "--plotFileFormat",
        "png"
    ])

    deeptools.plotFingerprint.main(args)

    return plotfile.name


def cleanup(*files):
    for file in files:
        if os.path.exists(file):
            os.remove(file)


def test_plotFingerprint_default():
    """Image comparison test for default plotFingerprint output."""

    args = (
        f"-b {TEST_DATA}test1.bam {TEST_DATA}test2.bam "
        "-l test1 test2"
    ).split()

    plotfile = run_plotFingerprint(args)

    try:
        res = compare_images(
            ROOT + "test_plotFingerprint_default.png",
            plotfile,
            tolerance
        )

        assert res is None, res

    finally:
        cleanup(plotfile)


def test_plotFingerprint_ggplot():
    """Image comparison test for --ggplot output."""

    args = (
        f"-b {TEST_DATA}test1.bam {TEST_DATA}test2.bam "
        "-l test1 test2 "
        "--ggplot"
    ).split()

    plotfile = run_plotFingerprint(args)

    try:
        res = compare_images(
            ROOT + "test_plotFingerprint_ggplot.png",
            plotfile,
            tolerance
        )

        assert res is None, res

    finally:
        cleanup(plotfile)


def test_plotFingerprint_quality_metrics_and_JSD():
    """
    Test --outQualityMetrics together with --JSDsample.
    """

    plotfile = NamedTemporaryFile(
        suffix=".png",
        prefix="deeptools_testfile_",
        delete=False
    )

    qcfile = NamedTemporaryFile(
        suffix=".tab",
        prefix="deeptools_testfile_",
        delete=False
    )

    args = (
        f"-b {TEST_DATA}test1.bam {TEST_DATA}test2.bam "
        f"-o {plotfile.name} "
        "--plotFileFormat png "
        "-l test1 test2 "
        f"--outQualityMetrics {qcfile.name} "
        f"--JSDsample {TEST_DATA}test1.bam"
    ).split()

    try:
        deeptools.plotFingerprint.main(args)

        with open(qcfile.name) as _foo:
            lines = [
                line.rstrip("\n").split("\t")
                for line in _foo
            ]

        assert len(lines) == 3, f"expected 3 lines, got {len(lines)}"

        header = lines[0]
        auc = header.index("AUC")
        jsd = header.index("JS Distance")

        rows = {
            row[0]: row
            for row in lines[1:]
        }

        assert abs(float(rows["test1"][auc]) - 0.39310288701202156) < 1e-4
        assert abs(float(rows["test2"][auc]) - 0.3641251150405128) < 1e-4
        assert abs(float(rows["test2"][jsd]) - 0.078613413909822) < 1e-4

    finally:
        cleanup(plotfile.name, qcfile.name)