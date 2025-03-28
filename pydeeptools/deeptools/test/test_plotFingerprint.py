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
    args = "-b {0}test1.bam {0}test2.bam -o {1} --plotFileFormat png".format(TEST_DATA, plotfile.name).split()
    deeptools.plotProfile.main(args)

    res = compare_images(ROOT + 'test_plotFingerprint_default.png', plotfile.name, tolerance)
    
    assert res is None, res
    
    os.remove(plotfile.name)