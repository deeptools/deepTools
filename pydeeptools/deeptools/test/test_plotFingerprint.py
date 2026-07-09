import os
import filecmp
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import deeptools.plotFingerprint

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_plotFingerprint/"

print (TEST_DATA)
print(ROOT)
tolerance = 13

def test_plotFingerprint_default():
    plotfile = NamedTemporaryFile(suffix='.png', prefix='deeptools_testfile_', delete=False)
    args = "-b {0}test1.bam {0}test2.bam -o {1} --plotFileFormat png -l test1 test2".format(TEST_DATA, plotfile.name).split()
    deeptools.plotFingerprint.main(args)

    res = compare_images(ROOT + 'test_plotFingerprint_default.png', plotfile.name, tolerance)

    assert res is None, res

    os.remove(plotfile.name)


def test_plotFingerprint_quality_metrics_and_JSD():
    """
    Test --outQualityMetrics together with --JSDsample. Neither the quality
    metrics table nor the JS-distance computation was covered elsewhere.
    """
    plotfile = NamedTemporaryFile(suffix='.png', prefix='deeptools_testfile_', delete=False)
    qcfile = NamedTemporaryFile(suffix='.tab', prefix='deeptools_testfile_', delete=False)
    args = ("-b {0}test1.bam {0}test2.bam -o {1} --plotFileFormat png -l test1 test2 "
            "--outQualityMetrics {2} --JSDsample {0}test1.bam".format(TEST_DATA, plotfile.name, qcfile.name)).split()
    deeptools.plotFingerprint.main(args)

    with open(qcfile.name) as _foo:
        lines = [line.rstrip("\n").split("\t") for line in _foo]

    # header + one row per sample
    assert len(lines) == 3, f"expected 3 lines, got {len(lines)}"
    header = lines[0]
    auc = header.index("AUC")
    jsd = header.index("JS Distance")

    rows = {row[0]: row for row in lines[1:]}
    assert abs(float(rows["test1"][auc]) - 0.39310288701202156) < 1e-4
    assert abs(float(rows["test2"][auc]) - 0.3641251150405128) < 1e-4
    # JS distance of the JSDsample (test1) against itself is nan; test2 is finite
    assert abs(float(rows["test2"][jsd]) - 0.078613413909822) < 1e-4

    os.remove(plotfile.name)
    os.remove(qcfile.name)