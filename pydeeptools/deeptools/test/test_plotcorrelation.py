import deeptools.plotCorrelation as pc

import os.path
from os import unlink
from matplotlib.testing.compare import compare_images
import tempfile
import pytest

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
COR_DATA_IN1 = ROOT + "multiBamSummary_result1.npz"
COR_PLOT_1 = ROOT + "plotCorrelation_result1.png"
COR_PLOT_2 = ROOT + "plotCorrelation_result2.png"


@pytest.mark.filterwarnings(
    "ignore:Attempting to set identical low and high xlims:UserWarning"
)
def test_correlation_plot_with_minimal_options():
    """
    Test minimal command line args for correlation plot
    with output as correlation matrix along with matrix
    """
    _, out_matrix = tempfile.mkstemp(suffix=".tsv")
    _, out_png = tempfile.mkstemp(suffix=".png")
    args = "--corData {} -p heatmap -c pearson -o {}  --outFileCorMatrix {}".format(COR_DATA_IN1, out_png, out_matrix).split()
    pc.main(args)

    _foo = open(out_matrix, "r")
    resp = _foo.readlines()[2]
    _foo.close()
    expected = "'bowtie2 test1.bam'\t1.0000\t1.0000\n"
    assert expected in resp, f"'{expected}' not found in '{resp}'"

    res = compare_images(COR_PLOT_1, out_png, 50)
    assert res is None, "Plots do not match"

def test_correlation_plot_scatter():
    """
    Test command line args for correlation plot with output as scatter plot
    """
    _, out_png2 = tempfile.mkstemp(suffix=".png")
    args = "--corData {} -p scatterplot -c pearson  -o {}".format(COR_DATA_IN1, out_png2).split()
    pc.main(args)

    res = compare_images(COR_PLOT_2, out_png2, 50)
    assert res is None, "Plots do not match"
