import os
import filecmp
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import deeptools.plotPCA

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_plotPCA/"

print(TEST_DATA)
print(ROOT)
tolerance = 50

def test_plotPCA_default():
    plotfile = NamedTemporaryFile(suffix='.png', prefix='deeptools_testfile_', delete=False)
    tsvfile  = NamedTemporaryFile(suffix='.tsv', prefix='deeptools_testfile_', delete=False)
    args = "-in {0}test_samples.npz -o {1} --outFileNameData {2}".format(TEST_DATA, plotfile.name, tsvfile.name).split()
    deeptools.plotPCA.main(args)

    res = compare_images(ROOT + 'test_plotPCA_default.png', plotfile.name, tolerance)
    assert res is None, res
    #assert filecmp.cmp(os.path.join(ROOT, 'test_plotPCA_default.tsv'), tsvfile.name) is True
    
    os.remove(plotfile.name)
    os.remove(tsvfile.name)