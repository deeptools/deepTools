import os
from tempfile import NamedTemporaryFile

import numpy as np
from matplotlib.testing.compare import compare_images

import deeptools.plotFingerprint

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_plotFingerprint/"
tolerance = 13


def run_plotFingerprint(args):
    """Run plotFingerprint and return generated plot file."""

    with NamedTemporaryFile(
        suffix=".png",
        prefix="deeptools_testfile_",
        delete=False
    ) as plotfile:
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

    with (
        NamedTemporaryFile(
            suffix=".png",
            prefix="deeptools_testfile_",
            delete=False
        ) as plotfile,
        NamedTemporaryFile(
            suffix=".tab",
            prefix="deeptools_testfile_",
            delete=False
        ) as qcfile,
    ):
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


def run_plotFingerprint_extendReads(bamfiles, labels, rawfile, qcfile):
    args = (
        ["-b"] + [TEST_DATA + x for x in bamfiles]
        + ["-l"] + labels
        + ["--extendReads",
           "--region", "chr2:4999000:5003000", "--binSize", "10", "--numberOfSamples", "400", "-p", "1",
           "--outRawCounts", rawfile, "--outQualityMetrics", qcfile]
    )
    deeptools.plotFingerprint.main(args)

    counts = np.loadtxt(rawfile, skiprows=2, ndmin=2)
    raw = {label: counts[:, idx] for idx, label in enumerate(labels)}
    with open(qcfile) as _foo:
        lines = [line.rstrip("\n").split("\t") for line in _foo]
    qc = {row[0]: np.array(row[1:], dtype=float) for row in lines[1:]}

    return raw, qc


def test_plotFingerprint_extendReads_single_bam(tmp_path):
    raw, qc = run_plotFingerprint_extendReads(["test_paired2.bam"], ["test_paired2"],
                                              str(tmp_path / "raw.txt"), str(tmp_path / "qc.txt"))
    expected_counts = np.loadtxt(ROOT + "test_plotFingerprint_extendReads_raw.txt", skiprows=2)
    np.testing.assert_array_equal(raw["test_paired2"], expected_counts)
    expected_qc = [0.10334787097692963, 0.4761528725688839, 0.595, 2.0396151481800195e-14,
                   0.745, 0.5373923764988261, 0.324362640621008]
    np.testing.assert_allclose(qc["test_paired2"], expected_qc, rtol=1e-6, atol=1e-12)


def test_plotFingerprint_extendReads_order_independent(tmp_path):
    bamfiles = ["test_paired.bam", "test_paired2.bam", "test_paired2.cram"]
    labels = ["paired", "paired2", "paired2_cram"]
    raw1, qc1 = run_plotFingerprint_extendReads(bamfiles, labels,
                                                str(tmp_path / "raw1.txt"), str(tmp_path / "qc1.txt"))
    raw2, qc2 = run_plotFingerprint_extendReads(bamfiles[::-1], labels[::-1],
                                                str(tmp_path / "raw2.txt"), str(tmp_path / "qc2.txt"))

    for label in labels:
        np.testing.assert_array_equal(raw1[label], raw2[label])
        np.testing.assert_array_equal(qc1[label], qc2[label])
