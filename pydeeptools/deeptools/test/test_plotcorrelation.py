import deeptools.plotCorrelation as pc

import os.path
from os import unlink


ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
COR_DATA_IN1 = ROOT + "multiBamSummary_result1.npz"
COR_PLOT_1 = ROOT + "plotCorrelation_result1.png"
COR_PLOT_2 = ROOT + "plotCorrelation_result2.png"


def test_correlation_plot_with_minimal_options():
    """
    Test minimal command line args for correlation plot
    with output as correlation matrix along with matrix
    """
    out_matrix = '/tmp/correlation_matrix.tsv'
    out_png = '/tmp/correlation_plot1.png'
    args = "--corData {} -p heatmap -c pearson -o {}  --outFileCorMatrix {}".format(COR_DATA_IN1, out_png, out_matrix).split()
    pc.main(args)

    _foo = open(out_matrix, "r")
    resp = _foo.readlines()[2]
    _foo.close()
    expected = "'bowtie2 test1.bam'\t1.0000\t1.0000\n"
    assert expected in resp, f"'{expected}' not found in '{resp}'"

    png_file_size = os.path.getsize(COR_PLOT_1)
    expected_file_size = os.path.getsize(out_png)
    size_tolerance = 5000
    size_difference = abs(png_file_size - expected_file_size)
    assert size_difference <= size_tolerance, "File size do not match"
    unlink(out_png)


def test_correlation_plot_scatter():
    """
    Test command line args for correlation plot with output as scatter plot
    """
    out_png2 = '/tmp/correlation_plot2.png'
    args = "--corData {} -p scatterplot -c pearson  -o {}".format(COR_DATA_IN1, out_png2).split()
    pc.main(args)
    png_file_size = os.path.getsize(COR_PLOT_2)
    expected_file_size = os.path.getsize(out_png2)
    size_tolerance = 5000
    size_difference = abs(png_file_size - expected_file_size)
    assert size_difference <= size_tolerance, "File size do not match"
    unlink(out_png2)
