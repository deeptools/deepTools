import os

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../../bin/"


def test_plot_coverage():
    os.system("{}plotCoverage --version".format(ROOT))


def test_bam_coverage():
    os.system("{}bamCoverage --version".format(ROOT))
