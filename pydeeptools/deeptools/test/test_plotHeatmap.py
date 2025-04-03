import os
import filecmp
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import deeptools.plotHeatmap
import deeptools.heatmapper_utilities

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_plotHeatmap/"

print (TEST_DATA)
print(ROOT)
tolerance = 13

def test_plotHeatmap_default():
    plotfile = NamedTemporaryFile(suffix='.png', prefix='deeptools_testfile_', delete=False)

    args = "--matrixFile {0}computeMatrix_result1.gz   --outFileName {1} ".format(TEST_DATA, plotfile.name).split()
    
    deeptools.plotHeatmap.main(args)

    # assert filecmp.cmp(os.path.join(ROOT, 'outMatrix_default.tab'), txtfile.name) is True
    res = compare_images(ROOT + 'plotHeatmap_default.png', plotfile.name, tolerance)
    assert res is None, res
    
    os.remove(plotfile.name)

