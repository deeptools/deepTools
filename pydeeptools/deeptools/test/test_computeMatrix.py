import math
import deeptools.computeMatrix2 as cm
import gzip
import tempfile
import os.path
import json
import numpy as np
from typing import Dict, List, Any, Tuple


ALLOWED_DELTA = 1.0
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
REGIONS_IN1 = ROOT + "computeMatrix1.bed"
REGIONS_IN2 = ROOT + "computeMatrix2.bed"
BIGWIG_IN1 = ROOT + "bamCoverage_result4.bw"
BIGWIG_IN2 = ROOT + "computeMatrix2.bw"
OUT_ARCHIEVE1 = ROOT + "computeMatrix_result1.gz"
OUT_ARCHIEVE2 = ROOT + "computeMatrix_result2.gz"
OUT_ARCHIEVE3 = ROOT + "computeMatrix_result3.gz"

def _parse_mat_gz(file_path: str) -> Tuple[Dict[str, Any], List[List[float]]]:
    """
    Parse a .mat.gz file into header and numeric data.

    Args:
        file_path: Path to the .mat.gz file

    Returns:
        Tuple of (header_dict, data_matrix) where data_matrix is list of lists of floats
    """
    with gzip.open(file_path, 'rt') as f:
        lines = f.readlines()

    if not lines:
        raise ValueError("File is empty")

    # Parse header (first line)
    header_line = lines[0].strip()
    if not header_line.startswith('@'):
        raise ValueError(f"Header must start with '@', got: {header_line[:20]}...")

    header_dict = json.loads(header_line[1:])

    # Parse data lines (starting from second line)
    data_lines = lines[1:]
    data_matrix = []

    for i, line in enumerate(data_lines):
        line = line.strip()
        if not line:
            continue

        parts = line.split('\t')
        # Skip header row if present
        if parts[0] == "chromosome":
            continue

        # Extract numeric values (from column 6 onwards)
        try:
            values = [float(x) if x != 'nan' else float('nan') for x in parts[6:]]
            data_matrix.append(values)
        except ValueError as e:
            raise ValueError(f"Error parsing numeric values in line {i+1}: {e}")

    return header_dict, data_matrix

def _compare_mat_gz(observed_file: str, expected_file: str) -> Tuple[List[str], List[str], Dict[str, List[str]]]:
    """
    Compare two .mat.gz files for header and data equality.
    """

    # Parse both files
    obs_header, obs_data = _parse_mat_gz(observed_file)
    exp_header, exp_data = _parse_mat_gz(expected_file)

    header_differences = []

    # Get all keys from both headers
    all_keys = set(obs_header.keys()) | set(exp_header.keys())

    for key in all_keys:
        if key not in obs_header:
            header_differences.append(f"Missing key in observed: {key}")
        elif key not in exp_header:
            header_differences.append(f"Missing key in expected: {key}")
        else:
            # Compare values
            obs_val = obs_header[key]
            exp_val = exp_header[key]

            # Special handling for NaN values
            if isinstance(obs_val, list) and isinstance(exp_val, list):
                # For lists, check if they're equal (handling NaN properly)
                if len(obs_val) != len(exp_val):
                    header_differences.append(f"List length mismatch for key '{key}': {len(obs_val)} vs {len(exp_val)}")
                else:
                    for i, (obs_item, exp_item) in enumerate(zip(obs_val, exp_val)):
                        if isinstance(obs_item, float) and isinstance(exp_item, float):
                            if np.isnan(obs_item) and np.isnan(exp_item):
                                continue  # Both NaN, equal
                            elif np.isnan(obs_item) or np.isnan(exp_item):
                                header_differences.append(f"NaN mismatch at index {i} in key '{key}': {obs_item} vs {exp_item}")
                            elif abs(obs_item - exp_item) > ALLOWED_DELTA:
                                header_differences.append(f"Float value mismatch at index {i} in key '{key}': {obs_item} vs {exp_item}")
                        elif obs_item != exp_item:
                            # Escape for group_labels, which get _r1, ... in python
                            if key == 'group_labels':
                                # if only 1 region, python implement defaults to 'genes'
                                if len(obs_val) == 1 and len(exp_val) == 1 and exp_val[0] == 'genes':
                                    continue
                                _obs_item = obs_item.split('/')[-1].split('.')[0]
                                _exp_item = exp_item.split('/')[-1].split('.')[0]
                                if _obs_item != _exp_item:
                                    header_differences.append(f"Value mismatch at index {i} in key '{key}': {_obs_item} vs {_exp_item}")
                            else:
                                header_differences.append(f"Value mismatch at index {i} in key '{key}': {obs_item} vs {exp_item}")
            else:
                # For non-list values, simple comparison
                if isinstance(obs_val, float) and isinstance(exp_val, float):
                    if np.isnan(obs_val) and np.isnan(exp_val):
                        continue  # Both NaN, equal
                    elif np.isnan(obs_val) or np.isnan(exp_val):
                        header_differences.append(f"NaN mismatch for key '{key}': {obs_val} vs {exp_val}")
                    elif abs(obs_val - exp_val) > ALLOWED_DELTA:
                        header_differences.append(f"Float value mismatch for key '{key}': {obs_val} vs {exp_val}")
                elif obs_val != exp_val:
                    header_differences.append(f"Value mismatch for key '{key}': {obs_val} vs {exp_val}")

    data_differences = []
    rowdiffdict = {}
    if len(obs_data) != len(exp_data):
        data_differences.append(f"Data rows mismatch: {len(obs_data)} vs {len(exp_data)}")
    else:
        # Check each row
        min_rows = min(len(obs_data), len(exp_data))
        for i in range(min_rows):
            obs_row = obs_data[i]
            exp_row = exp_data[i]

            if len(obs_row) != len(exp_row):
                data_differences.append(f"Row {i} column count mismatch: {len(obs_row)} vs {len(exp_row)}")
            else:
                # Check each value in the row
                for j in range(len(obs_row)):
                    obs_val = obs_row[j]
                    exp_val = exp_row[j]

                    if isinstance(obs_val, float) and isinstance(exp_val, float):
                        if np.isnan(obs_val) and np.isnan(exp_val):
                            continue  # Both NaN, equal
                        elif np.isnan(obs_val) or np.isnan(exp_val):
                            data_differences.append(f"NaN mismatch at row {i}, col {j}: {obs_val} vs {exp_val}")
                            if i not in rowdiffdict:
                                rowdiffdict[i] = [obs_row, exp_row]
                        elif abs(obs_val - exp_val) > ALLOWED_DELTA:
                            data_differences.append(f"Float value mismatch at row {i}, col {j}: {obs_val} vs {exp_val}")
                            if i not in rowdiffdict:
                                rowdiffdict[i] = [obs_row, exp_row]
                    elif obs_val != exp_val:
                        data_differences.append(f"Value mismatch at row {i}, col {j}: {obs_val} vs {exp_val}")
                        if i not in rowdiffdict:
                            rowdiffdict[i] = [obs_row, exp_row]

    return (header_differences, data_differences, rowdiffdict)

def _compare_tab_files(observed, expected, rtol=1e-3, atol=ALLOWED_DELTA):
    """Compare two tab files with floating-point tolerance."""
    with open(observed) as f: obs = [line.strip().split('\t') for line in f if line.strip() and not line.startswith('#')]
    with open(expected) as f: exp = [line.strip().split('\t') for line in f if line.strip() and not line.startswith('#')]

    if len(obs) != len(exp): return False, [f"Row count: {len(obs)} vs {len(exp)}"]

    diffs = []
    for i, (o, e) in enumerate(zip(obs, exp)):
        if len(o) != len(e):
            diffs.append(f"Row {i+1} col count: {len(o)} vs {len(e)}")
            continue
        for j, (a, b) in enumerate(zip(o, e)):
            try:
                fa, fb = float(a), float(b)

                if math.isnan(fa) or math.isnan(fb):
                    # Only equal if both are NaN
                    if not (math.isnan(fa) and math.isnan(fb)):
                        diffs.append(f"Row {i+1}, Col {j+1}: {a} vs {b}")
                elif not (fa == fb or (abs(fa-fb) <= atol + rtol*abs(fb))):
                    diffs.append(f"Row {i+1}, Col {j+1}: {a} vs {b}")
            except ValueError:
                # sames 'genes'/'regions' naming in python3 as before
                if ':' in a and ':' in b:
                    a = a.split(':')[1]
                    b = b.split(':')[1]
                # Deal with 'nan'
                if a != b:
                    diffs.append(f"Row {i+1}, Col {j+1}: '{a}' vs '{b}'")

    return len(diffs) == 0, diffs

def _compare_bed_files(observed, expected):
    """Compare BED files ignoring header lines and the last field."""
    with open(observed) as f:
        obs = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    with open(expected) as f:
        exp = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    diffs = []
    for i, (a, b) in enumerate(zip(obs, exp)):
        # Split lines and compare all fields except the last one
        a_fields = a.split('\t')
        b_fields = b.split('\t')

        if len(a_fields) != len(b_fields):
            diffs.append(f"Line {i+1}: field count mismatch {len(a_fields)} vs {len(b_fields)}")
        else:
            # Compare all fields except the last one
            if a_fields[:-1] != b_fields[:-1]:
                diffs.append(f"Line {i+1}: {a} vs {b} (ignoring last field)")
    if diffs:
        return False, diffs
    return True, []


def test_compute_matrix_refpoint():
    exp_npz = ROOT + "test_computematrix/mat.gz"
    exp_mat = ROOT + "test_computematrix/mat.tab"
    exp_bed = ROOT + "test_computematrix/mat.bed"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    _, outfile_bed = tempfile.mkstemp(suffix='.bed')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1}
    --scoreFileName {BIGWIG_IN1}
    -b 20 -a 20
    -o {outfile_npz} --outFileNameMatrix {outfile_mat} --outFileSortedRegions {outfile_bed}
    """.split()
    print(' '.join(args))
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

    # bed_match, bed_diffs = _compare_bed_files(outfile_bed, exp_bed)
    # assert bed_match, f"bed mismatch: {bed_diffs}"

def test_compute_matrix_refpoint_multibed():
    exp_npz = ROOT + "test_computematrix/mat_multibed.gz"
    exp_mat = ROOT + "test_computematrix/mat_multibed.tab"
    exp_bed = ROOT + "test_computematrix/mat_multibed.bed"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    _, outfile_bed = tempfile.mkstemp(suffix='.bed')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN1}
    --scoreFileName {BIGWIG_IN1}
    -b 20 -a 20
    -o {outfile_npz} --outFileNameMatrix {outfile_mat} --outFileSortedRegions {outfile_bed}
    """.split()
    print(' '.join(args))
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_multibw():
    exp_npz = ROOT + "test_computematrix/mat_multibw.gz"
    exp_mat = ROOT + "test_computematrix/mat_multibw.tab"
    exp_bed = ROOT + "test_computematrix/mat_multibw.bed"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    _, outfile_bed = tempfile.mkstemp(suffix='.bed')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN1}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN1} {BIGWIG_IN1} {BIGWIG_IN1}
    -b 20 -a 20
    -o {outfile_npz} --outFileNameMatrix {outfile_mat} --outFileSortedRegions {outfile_bed}
    """.split()
    print(' '.join(args))
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_a270b310():
    exp_npz = ROOT + "test_computematrix/mat_a270_b310.gz"
    exp_mat = ROOT + "test_computematrix/mat_a270_b310.tab"
    exp_bed = ROOT + "test_computematrix/mat_a270_b310.bed"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    _, outfile_bed = tempfile.mkstemp(suffix='.bed')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1}
    --scoreFileName {BIGWIG_IN1}
    -b 310 -a 270
    -o {outfile_npz} --outFileNameMatrix {outfile_mat} --outFileSortedRegions {outfile_bed}
    """.split()
    print(' '.join(args))
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

# def test_compute_matrix_with_reference_point_and_advance_options_1():
#     """
#     Test minimal command line args compute matrix with reference based
#     mode along with advance options of sorting using sum and average
#     type bin as sum
#     """
#     _, outfile = tempfile.mkstemp(suffix='.gz')
#     args = f"reference-point --regionsFileName {REGIONS_IN1} --scoreFileName {BIGWIG_IN1} -o {outfile} -bs 10 --sortUsing sum --averageTypeBins sum -b 10 -a 10".split()
#     print(' '.join(args))
#     cm.main(args)

#     archieve_file_size = os.path.getsize(OUT_ARCHIEVE1)
#     expected_file_size = os.path.getsize(outfile)
#     size_tolerance = 500
#     size_difference = abs(archieve_file_size - expected_file_size)
#     assert size_difference <= size_tolerance, "File size do not match"

# def test_compute_matrix_with_reference_point_and_advance_options_2():
#     """
#     Test minimal command line args compute matrix with reference based mode
#     with before and after region start length
#     """
#     _, outfile = tempfile.mkstemp(suffix='.gz')
#     args = "reference-point --regionsFileName {} --scoreFileName {} -o {} -bs 10 -b 10 -a 10".format(REGIONS_IN2, BIGWIG_IN2, outfile).split()
#     cm.main(args)

#     archieve_file_size = os.path.getsize(OUT_ARCHIEVE2)
#     expected_file_size = os.path.getsize(outfile)
#     size_tolerance = 500
#     size_difference = abs(archieve_file_size - expected_file_size)
#     assert size_difference <= size_tolerance, "File size do not match"

# def test_compute_matrix_with_scale_regions():
#     """
#     Test minimal command line args compute matrix with scale regions mode
#     """
#     _, outfile = tempfile.mkstemp(suffix='.gz')
#     args = "scale-regions --regionsFileName {} --scoreFileName {} -o {}".format(REGIONS_IN2, BIGWIG_IN2, outfile).split()
#     cm.main(args)

#     archieve_file_size = os.path.getsize(OUT_ARCHIEVE3)
#     expected_file_size = os.path.getsize(outfile)
#     size_tolerance = 500
#     size_difference = abs(archieve_file_size - expected_file_size)
#     assert size_difference <= size_tolerance, "File size do not match"
