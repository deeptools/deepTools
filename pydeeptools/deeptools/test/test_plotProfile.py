import os
import filecmp
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import deeptools.plotProfile
import deeptools.heatmapper_utilities

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_plotProfile/"

print (TEST_DATA)
print(ROOT)
tolerance = 13

def test_plotProfile_default():
    plotfile = NamedTemporaryFile(suffix='.png', prefix='deeptools_testfile_', delete=False)
    args = "--matrixFile {0}computeMatrix_result1.gz   --outFileName {1} --plotFileFormat png".format(TEST_DATA, plotfile.name).split()
    deeptools.plotProfile.main(args)

    res = compare_images(ROOT + 'plotProfile_default.png', plotfile.name, tolerance)
    
    assert res is None, res
    os.remove(plotfile.name)