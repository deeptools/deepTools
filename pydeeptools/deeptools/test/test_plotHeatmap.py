import os
import tempfile
import pytest
import matplotlib 
import matplotlib.pyplot as plt 
from matplotlib.testing.compare import compare_images
import deeptools.plotHeatmap

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_plotHeatmap/"

print(TEST_DATA)
print(ROOT)
tolerance = 13

@pytest.fixture(autouse=True)
def debug_matplotlib():
    print("\n--- matplotlib settings ---")
    print("backend:", matplotlib.get_backend())
    print("figure.figsize:", matplotlib.rcParams["figure.figsize"])
    print("figure.dpi:", matplotlib.rcParams["figure.dpi"])
    print("savefig.dpi:", matplotlib.rcParams["savefig.dpi"])
    print("font.size:", matplotlib.rcParams["font.size"])
    print("font.size:", matplotlib.rcParams["font.family"])
    yield

@pytest.fixture(autouse=True)
def cleanup_matplotlib():
    with matplotlib.rc_context():
        yield
    plt.close("all")


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
