"""Check exact row alignment across all computeMatrix output formats."""

import gzip
import json
import math

import numpy as np
import pyBigWig
import pytest

from deeptools import computeMatrix2


@pytest.fixture
def row_inputs(tmp_path):
    values = [[0, 0], [1, 3], [math.nan, 4], [7, -2],
              [math.nan, math.nan], [1, 3]]
    tracks = []
    for sample in range(2):
        path = tmp_path / f"sample{sample}.bw"
        with pyBigWig.open(str(path), "w") as bw:
            bw.addHeader([("chr1", 1000)])
            for index, row in enumerate(values):
                for column, value in enumerate(row):
                    if math.isfinite(value):
                        start = 100 + 100 * index + 5 * column
                        bw.addEntries(["chr1"], [start], ends=[start + 5],
                                      values=[float(value * (sample + 1))])
        tracks.append(str(path))
    groups = [[0], [2, 1, 5], [4, 3]]
    beds = []
    for index, group in enumerate(groups):
        bed = tmp_path / (["zeros", "signal", "missing"][index] + ".bed")
        bed.write_text("".join(
            f"chr1\t{100 + 100 * row}\t{110 + 100 * row}\trow{row}\t0\t+\n"
            for row in group))
        beds.append(str(bed))
    matrix = np.array([row + [2 * value for value in row] for row in values])
    return beds, tracks, groups, matrix


@pytest.mark.parametrize("mode", ["reference-point", "scale-regions"])
@pytest.mark.parametrize("sort,metric,samples", [
    ("keep", "mean", []), ("no", "mean", []),
    ("ascend", "mean", []), ("descend", "mean", []),
    ("ascend", "median", []), ("descend", "max", []),
    ("ascend", "min", []), ("descend", "sum", [2]),
    ("descend", "region_length", []),
])
@pytest.mark.parametrize("filtering", ["none", "zeros", "thresholds"])
def test_row_alignment(tmp_path, row_inputs, mode, sort, metric, samples, filtering):
    beds, tracks, groups, matrix = row_inputs
    gz = tmp_path / "matrix.gz"
    tab = tmp_path / "matrix.tab"
    bed = tmp_path / "sorted.bed"
    args = [mode, "-R", *beds, "-S", *tracks, "--binSize", "5", "-p", "1",
            "--sortRegions", sort, "--sortUsing", metric,
            "--samplesLabel", "one", "two", "-o", str(gz),
            "--outFileNameMatrix", str(tab), "--outFileSortedRegions", str(bed)]
    args += ["-b", "0", "-a", "10"] if mode == "reference-point" else ["-m", "10"]
    if samples:
        args += ["--sortUsingSamples", *map(str, samples)]
    if filtering == "zeros":
        args += ["--skipZeros"]
    elif filtering == "thresholds":
        args += ["--minThreshold", "0.5", "--maxThreshold", "10"]
    computeMatrix2.main(args)

    expected = []
    boundaries = [0]
    for group in groups:
        retained = []
        for index in group:
            row = matrix[index]
            if filtering == "zeros" and np.all(row == 0):
                continue
            if (filtering == "thresholds" and not np.isnan(row).any()
                    and (np.any(row <= 0.5) or np.any(row >= 10))):
                continue
            retained.append(index)
        if sort == "no":
            retained.sort()
        elif sort in ("ascend", "descend"):
            def key(index):
                if metric == "region_length":
                    return (False, 10)
                values = matrix[index, 2:] if samples else matrix[index]
                values = values[np.isfinite(values)]
                if len(values) == 0:
                    return (True, 0)
                return (False, getattr(np, metric)(values))
            retained.sort(key=key)
            if sort == "descend":
                retained.reverse()
        expected.extend(retained)
        boundaries.append(len(expected))

    with gzip.open(gz, "rt") as source:
        header = json.loads(source.readline()[1:])
        rows = [line.rstrip().split("\t") for line in source]
    assert header["group_boundaries"] == boundaries
    assert header["sample_boundaries"] == [0, 2, 4]
    assert header["group_labels"] == ["zeros", "signal", "missing"]
    assert [row[3] for row in rows] == [f"row{i}" for i in expected]
    observed = np.array([[float(x) for x in row[6:]] for row in rows])
    np.testing.assert_array_equal(observed, matrix[expected])
    raw_lines = tab.read_text().splitlines()
    assert raw_lines[0] == "#" + "\t".join(
        f"{label}:{boundaries[i + 1] - boundaries[i]}"
        for i, label in enumerate(["zeros", "signal", "missing"]))
    raw = np.array([[float(x) for x in line.split("\t")] for line in raw_lines[3:]])
    np.testing.assert_array_equal(raw, matrix[expected])
    bed_rows = [line.split("\t") for line in bed.read_text().splitlines()[1:]]
    assert [row[3] for row in bed_rows] == [f"row{i}" for i in expected]
    labels = {index: label for group, label in zip(groups, ["zeros", "signal", "missing"])
              for index in group}
    assert [row[-1] for row in bed_rows] == [labels[i] for i in expected]
