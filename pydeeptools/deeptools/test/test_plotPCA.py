import os
import filecmp
import numpy as np
import pytest
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
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
    full plotting path runs; set it False to exercise only the numeric output
    (e.g. small --ntop values whose plotting path is separately broken)."""
    tsvfile = NamedTemporaryFile(suffix='.tsv', prefix='deeptools_testfile_', delete=False)
    args = "-in {0}test_samples.npz --outFileNameData {1}".format(
        TEST_DATA, tsvfile.name).split()
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


def _sign_fix(coords, component_axis=1):
    """PCA eigenvector signs are arbitrary: they flip across BLAS/platforms and
    between implementations (e.g. sklearn vs a scipy/SVD rewrite). The sign
    freedom is per principal component, so normalize each component's vector to
    have a positive largest-magnitude entry.

    ``component_axis`` says which axis indexes the principal components:
    in the untransposed --outFileNameData table components are the columns
    (axis=1); in the transposed table they are the rows (axis=0)."""
    coords = np.array(coords, dtype=float)
    if component_axis == 1:
        for j in range(coords.shape[1]):
            i = np.argmax(np.abs(coords[:, j]))
            if coords[i, j] < 0:
                coords[:, j] = -coords[:, j]
    else:
        for i in range(coords.shape[0]):
            j = np.argmax(np.abs(coords[i, :]))
            if coords[i, j] < 0:
                coords[i, :] = -coords[i, :]
    return coords


# Golden eigenvalues captured from the sklearn-backed implementation over
# test_samples.npz with the default --ntop 500. Eigenvalues are the portable
# invariant (stable across BLAS backends and across the scipy/SVD rewrite);
# untransposed per-feature coordinates are not (see test_plotPCA_default_eigenvalues).
_GOLDEN_DEFAULT_EIGENVALUES = np.array([
    5.807692278756, 0.074230288836, 0.048971777735,
    0.036809415525, 0.026706723301, 0.017613563943,
])


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


def test_plotPCA_default_eigenvalues():
    """Regression on the untransposed eigenvalues, the portable numeric
    invariant of this path.

    We deliberately do NOT assert the projected coordinates here. After PC1
    the eigenvalues are tiny and near-degenerate (~0.07, 0.05, 0.04, ...), so
    the corresponding eigenvectors are free to rotate within that subspace,
    and the top-``ntop`` row selection (np.argpartition) breaks variance ties
    differently across BLAS backends (Linux OpenBLAS vs macOS Accelerate).
    The resulting per-feature coordinates are therefore not reproducible
    across platforms/implementations. Coordinate-level regression is covered
    by test_plotPCA_transpose, whose components are well separated and stable.
    The default plot itself is still pinned by test_plotPCA_default (image
    comparison)."""
    data = _run_pca()
    np.testing.assert_array_equal(data[:, 0], np.arange(1, 7))
    np.testing.assert_allclose(data[:, -1], _GOLDEN_DEFAULT_EIGENVALUES, rtol=1e-5)


def test_plotPCA_variance_matches_eigenvalues():
    """The per-PC variance fraction shown on the axis labels / scree plot is the
    eigenvalue proportion. Pin that relationship so the rewrite keeps the two in
    sync (eigenvalues are monotonically non-increasing and normalize to 1)."""
    eig = _run_pca()[:, -1]
    assert np.all(np.diff(eig) <= 1e-9), "eigenvalues must be non-increasing"
    pvar = eig / eig.sum()
    np.testing.assert_allclose(pvar.sum(), 1.0, rtol=1e-9)
    # PC1 dominates on this synthetic wt/kd matrix.
    assert pvar[0] > 0.9


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
    """When --ntop is below the sample count the table is truncated to the
    number of retained components (rows); every row still holds one loading
    per sample, so the plot of PC1 vs PC2 is well-defined."""
    data = _run_pca(["--ntop", "2"])
    assert data.shape == (2, 8)
    np.testing.assert_array_equal(data[:, 0], np.arange(1, 3))
    # First component carries all the variance for the 2-feature case.
    np.testing.assert_allclose(data[0, -1], 12.0, rtol=1e-6)
    assert abs(data[1, -1]) < 1e-6


def test_plotPCA_requested_PC_beyond_ntop_errors_cleanly():
    """Asking for a component that the retained rows cannot provide must exit
    with a clear message rather than crash with an IndexError."""
    plotfile = NamedTemporaryFile(suffix='.png', prefix='deeptools_testfile_', delete=False)
    args = "-in {0}test_samples.npz -o {1} --ntop 2 --PCs 1 3".format(TEST_DATA, plotfile.name).split()
    try:
        with pytest.raises(SystemExit) as exc:
            deeptools.plotPCA.main(args)
        assert "principal component" in str(exc.value)
    finally:
        if os.path.exists(plotfile.name):
            os.remove(plotfile.name)


def test_plotPCA_default_table_holds_sample_loadings():
    """In the default (untransposed) layout the table row for component i is
    the loading of every sample on PC i, i.e. row i of V^T from the SVD of the
    standardised (ntop rows x samples) matrix: the rows are orthonormal and
    equal, up to the per-component sign, to a numpy PCA on the same rows.
    (The bin scores U*S, one row per selected bin, are neither.)"""
    data = _run_pca()
    W = data[:, 1:7]
    np.testing.assert_allclose(W @ W.T, np.eye(6), atol=1e-6)

    mat = np.load(TEST_DATA + "test_samples.npz")["matrix"].astype(float)
    m = mat[np.argpartition(mat.var(axis=1), -500)[-500:]]
    m = (m - m.mean(axis=0)) / m.std(axis=0)
    _, _, Vt = np.linalg.svd(m - m.mean(axis=0), full_matrices=False)
    # PC1 dominates and is well separated, so its loadings are stable
    # across BLAS backends; compare it sign-invariantly.
    np.testing.assert_allclose(_sign_fix(W[:1].T)[:, 0], _sign_fix(Vt[:1].T)[:, 0], rtol=1e-5, atol=1e-8)


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
    args = "-in {0}test_samples.npz -o {1}".format(TEST_DATA, plotfile.name).split() + extra
    try:
        with pytest.raises(SystemExit) as exc:
            deeptools.plotPCA.main(args)
        assert msg in str(exc.value)
    finally:
        if os.path.exists(plotfile.name):
            os.remove(plotfile.name)


def test_plotPCA_requires_an_output():
    with pytest.raises(SystemExit) as exc:
        deeptools.plotPCA.main("-in {0}test_samples.npz".format(TEST_DATA).split())
    assert "must be specified" in str(exc.value)


# Golden values for the transposed PCA (samples as observations, so each row
# of the table is a component's projection across the six samples). Captured
# after fixing the projection bug; stored raw (components are rows -> axis=0).
_GOLDEN_TRANSPOSE_COORDS = np.array([
    [8.096369192617, 27.65422672360, -1.598082844166, -15.48892072797, -18.49767188707, -0.1659204570140],
    [3.552671394141, -4.837722476763, 20.02992087542, 0.4882670876827, -7.713434681130, -11.51970219935],
    [10.09925342625, -9.161879672887, 1.351881375625, -11.06368229223, -0.2099739792611, 8.984401142503],
    [-11.75933735644, 2.772538322120, 7.637287613484, -8.588903203498, 5.490108421837, 4.448306202499],
    [4.468893041249, 2.035996403779, -1.674415783401, -5.444924755301, 9.786060422166, -9.171609328491],
    [3.400996688227e-15, 3.400996688227e-15, 3.400996688227e-15, 3.400996688227e-15, 3.400996688227e-15, 3.400996688227e-15],
])
_GOLDEN_TRANSPOSE_EIGENVALUES = np.array([
    282.9918757435, 125.9323562327, 78.18623219731,
    65.59902453902, 47.29051128753, 1.388013416800e-29,
])


def test_plotPCA_transpose():
    """--transpose runs (previously crashed) and projects each sample onto the
    PCs. Coordinates are compared sign-invariantly; the last component is a
    numerical-zero residual so we skip its unstable sign."""
    data = _run_pca(["--transpose"])
    assert data.shape == (6, 8)
    np.testing.assert_array_equal(data[:, 0], np.arange(1, 7))
    # Transposed table: components are rows -> sign-fix per row (axis=0).
    coords = _sign_fix(data[:, 1:7], component_axis=0)
    golden = _sign_fix(_GOLDEN_TRANSPOSE_COORDS, component_axis=0)
    # Compare the informative components; the final ~1e-15 residual row is noise.
    np.testing.assert_allclose(coords[:-1], golden[:-1], rtol=1e-4, atol=1e-6)
    np.testing.assert_allclose(data[:, -1], _GOLDEN_TRANSPOSE_EIGENVALUES, rtol=1e-4, atol=1e-6)
    # Transposed eigenvalues differ from the untransposed layout.
    assert not np.allclose(data[:, -1], _GOLDEN_DEFAULT_EIGENVALUES)


def test_plotPCA_log2_and_rowCenter_affect_output():
    """--log2 and --rowCenter now actually transform the data before the PCA,
    so each changes the result relative to the default."""
    default = _run_pca()
    log2 = _run_pca(["--log2"])
    rowcenter = _run_pca(["--rowCenter"])
    assert not np.allclose(default[:, -1], log2[:, -1]), "--log2 was a no-op"
    assert not np.allclose(default[:, -1], rowcenter[:, -1]), "--rowCenter was a no-op"


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
        "-in {0}test_samples.npz "
        "-o {1} "
        "--outFileNameData {2} "
        "--ggplot"
    ).format(
        TEST_DATA,
        plotfile.name,
        tsvfile.name
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