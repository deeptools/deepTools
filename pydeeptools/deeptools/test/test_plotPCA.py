import os
import filecmp
import numpy as np
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


def test_plotPCA_outFileNameData():
    """
    Verify the numeric --outFileNameData output. The eigenvector sign is
    arbitrary (and can flip across BLAS/platforms), so we assert on the
    sign-independent eigenvalue column and the table shape rather than the
    projected coordinates.
    """
    plotfile = NamedTemporaryFile(suffix='.png', prefix='deeptools_testfile_', delete=False)
    tsvfile = NamedTemporaryFile(suffix='.tsv', prefix='deeptools_testfile_', delete=False)
    args = "-in {0}test_samples.npz -o {1} --outFileNameData {2}".format(TEST_DATA, plotfile.name, tsvfile.name).split()
    deeptools.plotPCA.main(args)

    # Columns: Component, wt1, wt2, wt3, kd1, kd2, kd3, Eigenvalue
    data = np.loadtxt(tsvfile.name, skiprows=1)
    assert data.shape == (6, 8), f"unexpected shape {data.shape}"
    # Component index column
    np.testing.assert_array_equal(data[:, 0], np.arange(1, 7))
    eigenvalues = data[:, -1]
    expected_eigenvalues = np.array([
        5.807692278755936, 0.07423028883557825, 0.04897177773493676,
        0.03680941552538939, 0.026706723301448194, 0.017613563942900697,
    ])
    np.testing.assert_allclose(eigenvalues, expected_eigenvalues, rtol=1e-5)

    os.remove(plotfile.name)
    os.remove(tsvfile.name)