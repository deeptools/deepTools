import os
from tempfile import NamedTemporaryFile

import numpy as np
import pytest
from matplotlib.testing.compare import compare_images

import deeptools.plotPCA

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_plotPCA/"

print(TEST_DATA)
print(ROOT)
tolerance = 50


def _run_pca(extra=None, plot=True):
    """Run plotPCA over the shared test matrix and return the parsed
    --outFileNameData table (header stripped). ``extra`` is a list of extra
    CLI tokens. When ``plot`` is True a plot file is also requested so the
    full plotting path runs; set it False to exercise only the numeric
    output."""
    tsvfile = NamedTemporaryFile(suffix='.tsv', prefix='deeptools_testfile_', delete=False)
    args = f"-in {TEST_DATA}test_samples.npz --outFileNameData {tsvfile.name}".split()
    plotfile = None
    if plot:
        plotfile = NamedTemporaryFile(suffix='.png', prefix='deeptools_testfile_', delete=False)
        args += ["-o", plotfile.name]
    if extra:
        args += extra
    deeptools.plotPCA.main(args)
    data = np.loadtxt(tsvfile.name, skiprows=1)
    os.remove(tsvfile.name)
    if plotfile is not None:
        os.remove(plotfile.name)
    return data


def _sign_fix(coords):
    coords = np.array(coords, dtype=float)
    for i in range(coords.shape[0]):
        j = np.argmax(np.abs(coords[i, :]))
        if coords[i, j] < 0:
            coords[i, :] = -coords[i, :]
    return coords

_GOLDEN_DEFAULT_COORDS = np.array([
    [8.096369192617, 27.65422672360, -1.598082844166, -15.48892072797, -18.49767188707, -0.1659204570140],
    [3.552671394141, -4.837722476763, 20.02992087542, 0.4882670876827, -7.713434681130, -11.51970219935],
    [10.09925342625, -9.161879672887, 1.351881375625, -11.06368229223, -0.2099739792611, 8.984401142503],
    [-11.75933735644, 2.772538322120, 7.637287613484, -8.588903203498, 5.490108421837, 4.448306202499],
    [4.468893041249, 2.035996403779, -1.674415783401, -5.444924755301, 9.786060422166, -9.171609328491],
    [3.400996688227e-15, 3.400996688227e-15, 3.400996688227e-15, 3.400996688227e-15, 3.400996688227e-15, 3.400996688227e-15],
])
_GOLDEN_DEFAULT_EIGENVALUES = np.array([
    282.9918757435, 125.9323562327, 78.18623219731,
    65.59902453902, 47.29051128753, 1.388013416800e-29,
])

def test_plotPCA_default():
    plotfile = NamedTemporaryFile(suffix='.png', prefix='deeptools_testfile_', delete=False)
    tsvfile  = NamedTemporaryFile(suffix='.tsv', prefix='deeptools_testfile_', delete=False)
    args = f"-in {TEST_DATA}test_samples.npz -o {plotfile.name} --outFileNameData {tsvfile.name}".split()
    deeptools.plotPCA.main(args)

    res = compare_images(ROOT + 'test_plotPCA_default.png', plotfile.name, tolerance)
    assert res is None, res

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
    args = f"-in {TEST_DATA}test_samples.npz -o {plotfile.name} --outFileNameData {tsvfile.name}".split()
    deeptools.plotPCA.main(args)

    # Columns: Component, wt1, wt2, wt3, kd1, kd2, kd3, Eigenvalue
    data = np.loadtxt(tsvfile.name, skiprows=1)
    assert data.shape == (6, 8), f"unexpected shape {data.shape}"
    # Component index column
    np.testing.assert_array_equal(data[:, 0], np.arange(1, 7))
    eigenvalues = data[:, -1]
    np.testing.assert_allclose(eigenvalues, _GOLDEN_DEFAULT_EIGENVALUES, rtol=1e-4, atol=1e-6)

    os.remove(plotfile.name)
    os.remove(tsvfile.name)


def test_plotPCA_default_eigenvalues():
    """Regression on the default (scores) eigenvalues and per-sample
    coordinates. Components are well separated and stable on this synthetic
    wt/kd matrix, so both are safe to pin sign-invariantly; the last
    component is a numerical-zero residual so we skip its unstable sign."""
    data = _run_pca()
    assert data.shape == (6, 8)
    np.testing.assert_array_equal(data[:, 0], np.arange(1, 7))
    # Components are rows -> sign-fix per row (axis=0).
    coords = _sign_fix(data[:, 1:7])
    golden = _sign_fix(_GOLDEN_DEFAULT_COORDS)
    # Compare the informative components; the final ~1e-15 residual row is noise.
    np.testing.assert_allclose(coords[:-1], golden[:-1], rtol=1e-4, atol=1e-6)
    np.testing.assert_allclose(data[:, -1], _GOLDEN_DEFAULT_EIGENVALUES, rtol=1e-4, atol=1e-6)


def test_plotPCA_variance_matches_eigenvalues():
    """The per-PC variance fraction shown on the axis labels / scree plot is the
    eigenvalue proportion. Pin that relationship so the rewrite keeps the two in
    sync (eigenvalues are monotonically non-increasing and normalize to 1)."""
    eig = _run_pca()[:, -1]
    assert np.all(np.diff(eig) <= 1e-9), "eigenvalues must be non-increasing"
    pvar = eig / eig.sum()
    np.testing.assert_allclose(pvar.sum(), 1.0, rtol=1e-9)
    # PC1 explains the largest (though not overwhelming) share of the
    # per-sample variance on this synthetic wt/kd matrix.
    assert pvar[0] > 0.4


def test_plotPCA_ntop_zero_uses_all_rows():
    """--ntop 0 disables the top-variable-rows filter and therefore changes the
    result relative to the default --ntop 500 (the test matrix has >500 rows)."""
    default = _run_pca()
    allrows = _run_pca(["--ntop", "0"])
    assert allrows.shape == (6, 8)
    # Different row selection -> different eigenvalues.
    assert not np.allclose(default[:, -1], allrows[:, -1])
    # Eigenvalues still normalize and stay ordered.
    eig = allrows[:, -1]
    assert np.all(np.diff(eig) <= 1e-9)


def test_plotPCA_ntop_smaller_than_samples():
    """When --ntop selects fewer bins than there are samples, the number of
    components is capped at the bin count (k = min(n_samples, n_bins)), so
    the table is truncated to that many rows; all 6 samples still appear as
    columns."""
    data = _run_pca(["--ntop", "2"], plot=False)
    assert data.shape == (2, 8)
    np.testing.assert_array_equal(data[:, 0], np.arange(1, 3))
    eig = data[:, -1]
    assert np.all(np.diff(eig) <= 1e-9), "eigenvalues must be non-increasing"
    assert np.all(eig >= -1e-9), "eigenvalues must be non-negative"
    # With exactly 2 standardized (unit population-variance) bins, the 2
    # components capture all the sample variance: sum of eigenvalues (which
    # use the n-1 denominator) equals 2 * n/(n-1) for n=6 samples.
    np.testing.assert_allclose(eig.sum(), 2 * 6 / 5, rtol=1e-6)


def test_plotPCA_ntop_below_samples_plots_successfully():
    """Plotting with fewer retained bins than samples still works: the number
    of available components is simply capped below the sample count
    (previously this crashed with an IndexError at correlation.py's scatter
    loop, which assumed one component per sample)."""
    plotfile = NamedTemporaryFile(suffix='.png', prefix='deeptools_testfile_', delete=False)
    args = f"-in {TEST_DATA}test_samples.npz -o {plotfile.name} --ntop 2".split()
    try:
        deeptools.plotPCA.main(args)
        assert os.path.exists(plotfile.name) and os.path.getsize(plotfile.name) > 0
    finally:
        if os.path.exists(plotfile.name):
            os.remove(plotfile.name)


def test_plotPCA_PCs_selection_does_not_change_table():
    """--PCs only selects which components are drawn; the numeric table always
    contains every component, so it is independent of --PCs."""
    default = _run_pca()
    pcs13 = _run_pca(["--PCs", "1", "3"])
    np.testing.assert_allclose(default, pcs13, rtol=1e-9, atol=1e-12)


@pytest.mark.parametrize("extra, msg", [
    (["--PCs", "2", "2"], "different principal components"),
    (["--PCs", "0", "1"], "at least 1"),
    (["--ntop", "-1"], "must be >= 0"),
])
def test_plotPCA_invalid_arguments_exit(extra, msg):
    plotfile = NamedTemporaryFile(suffix='.png', prefix='deeptools_testfile_', delete=False)
    args = f"-in {TEST_DATA}test_samples.npz -o {plotfile.name}".split() + extra
    try:
        with pytest.raises(SystemExit) as exc:
            deeptools.plotPCA.main(args)
        assert msg in str(exc.value)
    finally:
        if os.path.exists(plotfile.name):
            os.remove(plotfile.name)


def test_plotPCA_requires_an_output():
    with pytest.raises(SystemExit) as exc:
        deeptools.plotPCA.main(f"-in {TEST_DATA}test_samples.npz".split())
    assert "must be specified" in str(exc.value)


def test_plotPCA_log2_affects_output():
    """--log2 transforms the data before the PCA, so it changes the result
    relative to the default."""
    default = _run_pca()
    log2 = _run_pca(["--log2"])
    assert not np.allclose(default[:, -1], log2[:, -1]), "--log2 was a no-op"


def test_plotPCA_ggplot():
    """Image comparison test for --ggplot output."""

    plotfile = NamedTemporaryFile(
        suffix='.png',
        prefix='deeptools_testfile_',
        delete=False
    )

    tsvfile = NamedTemporaryFile(
        suffix='.tsv',
        prefix='deeptools_testfile_',
        delete=False
    )

    args = (
        f"-in {TEST_DATA}test_samples.npz "
        f"-o {plotfile.name} "
        f"--outFileNameData {tsvfile.name} "
        "--ggplot"
    ).split()

    try:
        deeptools.plotPCA.main(args)

        res = compare_images(
            ROOT + "test_plotPCA_ggplot.png",
            plotfile.name,
            tolerance
        )

        assert res is None, res

    finally:
        os.remove(plotfile.name)
        os.remove(tsvfile.name)
