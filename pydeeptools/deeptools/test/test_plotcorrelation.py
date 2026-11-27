import os.path
import tempfile

import numpy as np
import pytest
from matplotlib.testing.compare import compare_images

import deeptools.plotCorrelation as pc
from deeptools.correlation import Correlation

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
COR_DATA_IN1 = ROOT + "multiBamSummary_result1.npz"
COR_PLOT_1 = ROOT + "plotCorrelation_result1.png"
COR_PLOT_2 = ROOT + "plotCorrelation_result2.png"

COR_PLOT_GG_1 = ROOT + "plotCorrelation_result1_ggplot.png"
COR_PLOT_GG_2 = ROOT + "plotCorrelation_result2_ggplot.png"


def test_correlation_labels_length_mismatch_exits():
    _, out_png = tempfile.mkstemp(suffix=".png")
    args = f"--corData {COR_DATA_IN1} -p heatmap -c pearson -o {out_png} " \
           "--labels sample1 sample2 sample3".split()
    with pytest.raises(SystemExit):
        pc.main(args)


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
    args = f"--corData {COR_DATA_IN1} -p heatmap -c pearson -o {out_png}  --outFileCorMatrix {out_matrix}".split()
    pc.main(args)

    with open(out_matrix, "r") as _foo:
        resp = _foo.readlines()[2]
    expected = "'bowtie2 test1.bam'\t1.0000\t1.0000\n"
    assert expected in resp, f"'{expected}' not found in '{resp}'"

    res = compare_images(COR_PLOT_1, out_png, 50)
    assert res is None, "Plots do not match"

def test_correlation_plot_scatter():
    """
    Test command line args for correlation plot with output as scatter plot
    """
    _, out_png2 = tempfile.mkstemp(suffix=".png")
    args = f"--corData {COR_DATA_IN1} -p scatterplot -c pearson  -o {out_png2}".split()
    pc.main(args)

    res = compare_images(COR_PLOT_2, out_png2, 50)
    assert res is None, "Plots do not match"

def test_correlation_plot_with_minimal_options_ggplot():
    """
    Test minimal command line args for correlation plot
    with output as correlation matrix along with matrix
    """
    _, out_matrix_gg = tempfile.mkstemp(suffix=".tsv")
    _, out_png_gg = tempfile.mkstemp(suffix=".png")
    args = f"--corData {COR_DATA_IN1} -p heatmap -c pearson -o {out_png_gg}  --outFileCorMatrix {out_matrix_gg} --ggplot".split()
    pc.main(args)

    with open(out_matrix_gg, "r") as _foo:
        resp = _foo.readlines()[2]
    expected = "'bowtie2 test1.bam'\t1.0000\t1.0000\n"
    assert expected in resp, f"'{expected}' not found in '{resp}'"

    res = compare_images(COR_PLOT_GG_1, out_png_gg, 50)
    assert res is None, "Plots do not match"

def test_correlation_plot_scatter_ggplot():
    """
    Test command line args for correlation plot with output as scatter plot
    """
    _, out_png2_gg = tempfile.mkstemp(suffix=".png")
    args = f"--corData {COR_DATA_IN1} -p scatterplot -c pearson  -o {out_png2_gg} --ggplot".split()
    pc.main(args)

    res = compare_images(COR_PLOT_GG_2, out_png2_gg, 50)
    assert res is None, "Plots do not match"


def test_get_outlier_indices_uses_deviations_from_median():
    # Baseline centered around 100 with unit spread plus one extreme outlier at 1000.
    # Under Iglewicz & Hoaglin (1993), MAD = median(|x - median(x)|) = 1.0,
    # giving modified z-score |1000 - 100| / (1.4826 * 1.0) ~= 607 > 200.
    base = np.array([99.0, 100.0, 101.0] * 30 + [1000.0])
    outliers = Correlation.get_outlier_indices(base)
    assert list(outliers) == [len(base) - 1]

    # Modified Z-score outlier detection must be invariant to additive shifts.
    shifted = base + 500.0
    assert list(Correlation.get_outlier_indices(shifted)) == list(outliers)

