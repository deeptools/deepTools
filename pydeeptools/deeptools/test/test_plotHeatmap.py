import os
import tempfile
from matplotlib.testing.compare import compare_images
import deeptools.plotHeatmap

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_plotHeatmap/"

print (TEST_DATA)
print(ROOT)
tolerance = 13

def test_plotHeatmap_default():
    _, plotfile = tempfile.mkstemp(suffix=".png")

    args = f"--matrixFile {TEST_DATA}computeMatrix_result1.gz --outFileName {plotfile}".split()

    deeptools.plotHeatmap.main(args)

    # assert filecmp.cmp(os.path.join(ROOT, 'outMatrix_default.tab'), txtfile.name) is True
    res = compare_images(ROOT + 'plotHeatmap_default.png', plotfile, tolerance)
    assert res is None, f"{plotfile} doesn't match {ROOT}plotHeatmap_default.png"
