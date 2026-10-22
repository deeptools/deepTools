import glob
import gzip
import os

import numpy as np
import pytest

from deeptools import heatmapper

ROOT = os.path.dirname(os.path.abspath(__file__))


def _matrix_files():
    files = []
    for f in sorted(glob.glob(os.path.join(ROOT, "**", "*.gz"), recursive=True)):
        with gzip.open(f) as fh:
            if fh.read(1) == b"@":
                files.append(f)
    return files


MATRIX_FILES = _matrix_files()


def _read_fast(path, monkeypatch):
    # make sure the Rust reader is used and does not fall back to the per-line reader
    def no_fallback(self, matrix_file):
        raise AssertionError(f"fell back to the per-line reader for {matrix_file}")

    hm = heatmapper.heatmapper()
    with monkeypatch.context() as m:
        m.setattr(heatmapper.heatmapper, "_read_matrix_file_legacy", no_fallback)
        hm.read_matrix_file(path)
    return hm


def _read_legacy(path):
    hm = heatmapper.heatmapper()
    hm._read_matrix_file_legacy(path)
    return hm


def assert_same_heatmapper(a, b):
    ma, mb = a.matrix, b.matrix
    assert type(ma.matrix) is type(mb.matrix)
    assert ma.matrix.dtype == mb.matrix.dtype
    assert ma.matrix.shape == mb.matrix.shape
    assert np.array_equal(np.ma.getmaskarray(ma.matrix), np.ma.getmaskarray(mb.matrix))
    # compare the bits, so that NaN and -0.0 are checked as well
    da = np.ascontiguousarray(np.ma.getdata(ma.matrix))
    db = np.ascontiguousarray(np.ma.getdata(mb.matrix))
    assert np.array_equal(da.view(np.int64), db.view(np.int64))
    assert ma.regions == mb.regions
    assert ma.group_boundaries == mb.group_boundaries
    assert ma.sample_boundaries == mb.sample_boundaries
    assert ma.group_labels == mb.group_labels
    assert ma.sample_labels == mb.sample_labels
    assert ma.sort_method == mb.sort_method
    assert ma.sort_using == mb.sort_using
    assert a.parameters == b.parameters


def test_found_matrix_files():
    assert len(MATRIX_FILES) > 50


@pytest.mark.parametrize("path", MATRIX_FILES, ids=lambda p: os.path.relpath(p, ROOT))
def test_read_matrix_file_matches_legacy(path, monkeypatch):
    assert_same_heatmapper(_read_fast(path, monkeypatch), _read_legacy(path))


def test_read_matrix_file_fallback(tmp_path):
    src = os.path.join(ROOT, "test_data", "computeMatrixOperations.mat.gz")
    expected = _read_legacy(src)

    # file objects are read by the per-line reader
    hm = heatmapper.heatmapper()
    with open(src, "rb") as fh:
        hm.read_matrix_file(fh)
    assert_same_heatmapper(hm, expected)

    # so are values that only Python's float() accepts ("1_0" == 10.0)
    with gzip.open(src, "rt") as fh:
        lines = fh.read().split("\n")
    fields = lines[1].split("\t")
    fields[6] = "1_0"
    lines[1] = "\t".join(fields)
    odd = tmp_path / "odd.mat.gz"
    with gzip.open(odd, "wt") as fh:
        fh.write("\n".join(lines))
    hm = heatmapper.heatmapper()
    hm.read_matrix_file(str(odd))
    assert hm.matrix.matrix[0, 0] == 10.0
    assert_same_heatmapper(hm, _read_legacy(str(odd)))
