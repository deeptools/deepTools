import math
import deeptools.computeMatrix2 as cm
import gzip
import tempfile
import os.path
import json
import numpy as np
from typing import Dict, List, Any, Tuple


ALLOWED_DELTA = 1.0
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/test_computematrix/"
REGIONS_IN1 = ROOT + "input_computeMatrix_regions1.bed"
REGIONS_IN2 = ROOT + "input_computeMatrix_regions2.bed"
REGIONS_GTF = ROOT + "input_computeMatrix_regions3.gtf"
REGIONS_BED12 = ROOT + "input_computeMatrix_regions4bed12.bed"
BL = ROOT + "input_computeMatrix_blacklist.bed"
BIGWIG_IN1 = ROOT + "input_computeMatrix_bw1.bw"
BIGWIG_IN2 = ROOT + "input_computeMatrix_bw2.bw"
BIGWIG_IN3 = ROOT + "input_computeMatrix_bw3.bw"
BIGWIG_IN4 = ROOT + "input_computeMatrix_bw4.bw"

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
            # scale-regions, no default has no 'startlabel' or 'endlabel' in json header
            if key == 'startLabel' or key == 'endLabel':
                continue
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
    exp_npz = ROOT + "mat.gz"
    exp_mat = ROOT + "mat.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1}
    --scoreFileName {BIGWIG_IN1}
    -b 20 -a 20
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint2():
    exp_npz = ROOT + "mat2.gz"
    exp_mat = ROOT + "mat2.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 30 -a 30 --binSize 30
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_multibed():
    exp_npz = ROOT + "mat_multibed.gz"
    exp_mat = ROOT + "mat_multibed.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1}
    -b 20 -a 20
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_multibw():
    exp_npz = ROOT + "mat_multibw.gz"
    exp_mat = ROOT + "mat_multibw.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 20 -a 20
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_a270b310():
    exp_npz = ROOT + "mat_a270_b310.gz"
    exp_mat = ROOT + "mat_a270_b310.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1}
    --scoreFileName {BIGWIG_IN1}
    -b 310 -a 270
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_bs30():
    exp_npz = ROOT + "mat_bs30.gz"
    exp_mat = ROOT + "mat_bs30.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1}
    --scoreFileName {BIGWIG_IN1}
    -b 210 -a 210 --binSize 30
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_bed():
    exp_npz = ROOT + "mat_bs30.gz"
    exp_bed = ROOT + "mat_bs30.bed"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_bed = tempfile.mkstemp(suffix='.bed')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1}
    --scoreFileName {BIGWIG_IN1}
    -b 210 -a 210 --binSize 30
    -o {outfile_npz} --outFileSortedRegions {outfile_bed}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    bed_match, bed_diffs = _compare_bed_files(outfile_bed, exp_bed)
    assert bed_match, f"bed mismatch: {bed_diffs}"

def test_compute_matrix_refpoint_TES():
    exp_npz = ROOT + "mat_tes.gz"
    exp_mat = ROOT + "mat_tes.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1}
    --scoreFileName {BIGWIG_IN1}
    -b 210 -a 210 --binSize 30
    --referencePoint TES
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_center():
    exp_npz = ROOT + "mat_center.gz"
    exp_mat = ROOT + "mat_center.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1}
    --scoreFileName {BIGWIG_IN1}
    -b 210 -a 210 --binSize 30
    --referencePoint center
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_nanafterend():
    exp_npz = ROOT + "mat_nan.gz"
    exp_mat = ROOT + "mat_nan.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1}
    --scoreFileName {BIGWIG_IN1}
    -b 210 -a 210 --binSize 30
    --nanAfterEnd
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_sort_no():
    exp_npz = ROOT + "mat_sort_no.gz"
    exp_mat = ROOT + "mat_sort_no.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1}
    --scoreFileName {BIGWIG_IN1}
    -b 210 -a 210 --binSize 30
    --sortRegions no
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_sort_ascend():
    exp_npz = ROOT + "mat_ascend.gz"
    exp_mat = ROOT + "mat_ascend.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1}
    --scoreFileName {BIGWIG_IN1}
    -b 210 -a 210 --binSize 30
    --sortRegions ascend
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_sort_descend():
    exp_npz = ROOT + "mat_descend.gz"
    exp_mat = ROOT + "mat_descend.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1}
    --scoreFileName {BIGWIG_IN1}
    -b 210 -a 210 --binSize 30
    --sortRegions descend
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_sort_descend_median():
    exp_npz = ROOT + "mat_descend_median.gz"
    exp_mat = ROOT + "mat_descend_median.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 210 -a 210 --binSize 30
    --sortRegions descend
    --sortUsing median
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_sort_descend_max():
    exp_npz = ROOT + "mat_descend_max.gz"
    exp_mat = ROOT + "mat_descend_max.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 210 -a 210 --binSize 30
    --sortRegions descend
    --sortUsing max
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_sort_descend_min():
    exp_npz = ROOT + "mat_descend_min.gz"
    exp_mat = ROOT + "mat_descend_min.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 210 -a 210 --binSize 30
    --sortRegions descend
    --sortUsing min
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_sort_descend_sum():
    exp_npz = ROOT + "mat_descend_sum.gz"
    exp_mat = ROOT + "mat_descend_sum.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 210 -a 210 --binSize 30
    --sortRegions descend
    --sortUsing sum
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_sort_descend_regionlength():
    exp_npz = ROOT + "mat_descend_regionlength.gz"
    exp_mat = ROOT + "mat_descend_regionlength.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 210 -a 210 --binSize 30
    --sortRegions descend
    --sortUsing region_length
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_sort_usingsamples():
    exp_npz = ROOT + "mat_sortusingsamples.gz"
    exp_mat = ROOT + "mat_sortusingsamples.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 210 -a 210 --binSize 30
    --sortRegions descend
    --sortUsing sum
    --sortUsingSamples 1 3
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_avgtypemedian():
    exp_npz = ROOT + "mat_avgtypemedian.gz"
    exp_mat = ROOT + "mat_avgtypemedian.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 210 -a 210 --binSize 30
    --averageTypeBins median
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_avgtypemin():
    exp_npz = ROOT + "mat_avgtypemin.gz"
    exp_mat = ROOT + "mat_avgtypemin.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 210 -a 210 --binSize 30
    --averageTypeBins min
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_avgtypemax():
    exp_npz = ROOT + "mat_avgtypemax.gz"
    exp_mat = ROOT + "mat_avgtypemax.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 210 -a 210 --binSize 30
    --averageTypeBins max
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_avgtypestd():
    exp_npz = ROOT + "mat_avgtypestd.gz"
    exp_mat = ROOT + "mat_avgtypestd.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 210 -a 210 --binSize 30
    --averageTypeBins std
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_avgtypesum():
    exp_npz = ROOT + "mat_avgtypesum.gz"
    exp_mat = ROOT + "mat_avgtypesum.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 210 -a 210 --binSize 30
    --averageTypeBins sum
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_missingdataaszero():
    exp_npz = ROOT + "mat_maz.gz"
    exp_mat = ROOT + "mat_maz.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 30 -a 30 --binSize 30
    --missingDataAsZero
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_skipzeros():
    exp_npz = ROOT + "mat_skipzeros.gz"
    exp_mat = ROOT + "mat_skipzeros.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 210 -a 210 --binSize 30
    --skipZeros
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_thresh():
    exp_npz = ROOT + "mat_thresh.gz"
    exp_mat = ROOT + "mat_thresh.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 210 -a 210 --binSize 30
    --minThreshold 10000 --maxThreshold 50000
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_bl():
    exp_npz = ROOT + "mat_bl.gz"
    exp_mat = ROOT + "mat_bl.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 30 -a 30 --binSize 30
    --blackListFileName {BL}
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_sampleslabel():
    exp_npz = ROOT + "mat_sampleslabel.gz"
    exp_mat = ROOT + "mat_sampleslabel.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 30 -a 30 --binSize 30
    --samplesLabel "a" "b" "c" "d"
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_scale():
    exp_npz = ROOT + "mat_scale.gz"
    exp_mat = ROOT + "mat_scale.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 30 -a 30 --binSize 30
    --scale 3
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_bed12():
    exp_npz = ROOT + "mat_bed12.gz"
    exp_mat = ROOT + "mat_bed12.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_BED12}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 30 -a 30 --binSize 30
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_gtf():
    exp_npz = ROOT + "mat_gtf.gz"
    exp_mat = ROOT + "mat_gtf.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 30 -a 30 --binSize 30
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_bed12mg():
    exp_npz = ROOT + "mat_bed12mg.gz"
    exp_mat = ROOT + "mat_bed12mg.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_BED12}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 30 -a 30 --binSize 30
    --metagene
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_gtfmg():
    exp_npz = ROOT + "mat_gtfmg.gz"
    exp_mat = ROOT + "mat_gtfmg.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 30 -a 30 --binSize 30
    --metagene
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_gtfmglarge():
    exp_npz = ROOT + "mat_gtfmglarge.gz"
    exp_mat = ROOT + "mat_gtfmglarge.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 2000 -a 2000 --binSize 10
    --metagene
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_gtflarge():
    exp_npz = ROOT + "mat_gtflarge.gz"
    exp_mat = ROOT + "mat_gtflarge.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 2000 -a 2000 --binSize 10
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_refpoint_gtfmgid():
    exp_npz = ROOT + "mat_gtfmg_id.gz"
    exp_mat = ROOT + "mat_gtfmg_id.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    reference-point
    --regionsFileName {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3}
    -b 30 -a 30 --binSize 30
    --metagene
    --transcriptID gene
    --exonID CDS
    --transcript_id_designator gene_id
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"
