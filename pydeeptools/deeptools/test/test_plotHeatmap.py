import os
import tempfile
import matplotlib.pyplot as plt
from matplotlib.testing.compare import compare_images
import deeptools.plotHeatmap

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_plotHeatmap/"

tolerance = 13


def test_plotHeatmap_default():

    _, plotfile = tempfile.mkstemp(suffix=".png")

    args = (
        f"--matrixFile {TEST_DATA}computeMatrix_result1.gz "
        f"--outFileName {plotfile}"
    ).split()

    # Turn off interactive plotting
    plt.ioff()

    deeptools.plotHeatmap.main(args)

    res = compare_images(ROOT + "plotHeatmap_default.png", plotfile, tolerance)
    

    # Clean up
    try:
        os.remove(plotfile)
    except OSError:
        pass

    assert res is None, f"{plotfile} doesn't match {ROOT}plotHeatmap_default.png"
