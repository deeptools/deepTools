import deeptools.computeMatrix2 as cm
import tempfile
import os.path
import json
import numpy as np
from typing import Dict, List, Any, Tuple
from .test_computeMatrix_referencepoint import _parse_mat_gz, _compare_mat_gz, _compare_tab_files,_compare_bed_files

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

def test_compute_matrix_scaleregions():
    exp_npz = ROOT + "srmat.gz"
    exp_mat = ROOT + "srmat.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
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

def test_compute_matrix_scaleregions_multireg():
    exp_npz = ROOT + "srmat_multireg.gz"
    exp_mat = ROOT + "srmat_multireg.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 20 -a 20
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_multireg2():
    exp_npz = ROOT + "srmat_multireg2.gz"
    exp_mat = ROOT + "srmat_multireg2.tab"
    exp_bed = ROOT + "srmat_multireg2.bed"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    _, outfile_bed = tempfile.mkstemp(suffix='.bed')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 150 -a 195 --binSize 15 -m 150
    -o {outfile_npz} --outFileNameMatrix {outfile_mat} --outFileSortedRegions {outfile_bed}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

    bed_match, bed_diffs = _compare_bed_files(outfile_bed, exp_bed)
    assert bed_match, f"bed mismatch: {bed_diffs}"

def test_compute_matrix_scaleregions_reglabels():
    exp_npz = ROOT + "srmat_reglabel.gz"
    exp_mat = ROOT + "srmat_reglabel.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 500 -a 800 --binSize 20 -m 2400
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    --startLabel fakepeak_start
    --endLabel fakepeak_end
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_unscaled():
    exp_npz = ROOT + "srmat_unscaled.gz"
    exp_mat = ROOT + "srmat_unscaled.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 500 -a 800 --binSize 20 -m 2400
    --unscaled5prime 240 --unscaled3prime 40
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_unscaled2():
    exp_npz = ROOT + "srmat_unscaled2.gz"
    exp_mat = ROOT + "srmat_unscaled2.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 500 -a 800 --binSize 20 -m 200
    --unscaled5prime 20 --unscaled3prime 120
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_unscaled3():
    exp_npz = ROOT + "srmat_unscaled3.gz"
    exp_mat = ROOT + "srmat_unscaled3.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"


def test_compute_matrix_scaleregions_median():
    exp_npz = ROOT + "srmat_med.gz"
    exp_mat = ROOT + "srmat_med.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --averageTypeBins median
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_min():
    exp_npz = ROOT + "srmat_min.gz"
    exp_mat = ROOT + "srmat_min.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --averageTypeBins min
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_std():
    exp_npz = ROOT + "srmat_std.gz"
    exp_mat = ROOT + "srmat_std.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --averageTypeBins std
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_max():
    exp_npz = ROOT + "srmat_max.gz"
    exp_mat = ROOT + "srmat_max.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --averageTypeBins max
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"


def test_compute_matrix_scaleregions_sum():
    exp_npz = ROOT + "srmat_sum.gz"
    exp_mat = ROOT + "srmat_sum.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --averageTypeBins sum
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_descend():
    exp_npz = ROOT + "srmat_desc.gz"
    exp_mat = ROOT + "srmat_desc.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --sortRegions descend
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_ascend():
    exp_npz = ROOT + "srmat_ascend.gz"
    exp_mat = ROOT + "srmat_ascend.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --sortRegions ascend
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_sortno():
    exp_npz = ROOT + "srmat_no.gz"
    exp_mat = ROOT + "srmat_no.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --sortRegions no
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_sortusum():
    exp_npz = ROOT + "srmat_sortusum.gz"
    exp_mat = ROOT + "srmat_sortusum.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --sortRegions descend
    --sortUsing sum
    --sortUsingSamples 3
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_sortumin():
    exp_npz = ROOT + "srmat_sortumin.gz"
    exp_mat = ROOT + "srmat_sortumin.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --sortRegions descend
    --sortUsing min
    --sortUsingSamples 3
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_sortumedian():
    exp_npz = ROOT + "srmat_sortumedian.gz"
    exp_mat = ROOT + "srmat_sortumedian.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --sortRegions descend
    --sortUsing median
    --sortUsingSamples 3
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_sortumean():
    exp_npz = ROOT + "srmat_sortumean.gz"
    exp_mat = ROOT + "srmat_sortumean.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --sortRegions descend
    --sortUsing mean
    --sortUsingSamples 3
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_sortumax():
    exp_npz = ROOT + "srmat_sortumax.gz"
    exp_mat = ROOT + "srmat_sortumax.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --sortRegions descend
    --sortUsing max
    --sortUsingSamples 3
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_sortulen():
    exp_npz = ROOT + "srmat_sortulen.gz"
    exp_mat = ROOT + "srmat_sortulen.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --sortRegions descend
    --sortUsing region_length
    --sortUsingSamples 3
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_missingdataaszero():
    exp_npz = ROOT + "srmat_missingaszero.gz"
    exp_mat = ROOT + "srmat_missingaszero.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --sortRegions descend
    --sortUsing region_length
    --sortUsingSamples 3
    --missingDataAsZero
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_skipzeros():
    exp_npz = ROOT + "srmat_skipzeros.gz"
    exp_mat = ROOT + "srmat_skipzeros.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --sortRegions descend
    --sortUsing region_length
    --sortUsingSamples 3
    --missingDataAsZero
    --skipZeros
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_thresh():
    exp_npz = ROOT + "srmat_thresh.gz"
    exp_mat = ROOT + "srmat_thresh.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --sortRegions descend
    --sortUsing region_length
    --sortUsingSamples 3
    --minThreshold 12000 --maxThreshold 25000
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_bl():
    exp_npz = ROOT + "srmat_bl.gz"
    exp_mat = ROOT + "srmat_bl.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 180 -a 180 --binSize 20 -m 1000
    --unscaled5prime 20 --unscaled3prime 20
    --sortRegions descend
    --sortUsing region_length
    --blackListFileName {BL}
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_label():
    exp_npz = ROOT + "srmat_label.gz"
    exp_mat = ROOT + "srmat_label.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 70 -a 140 --binSize 35 -m 420
    --unscaled5prime 35 --unscaled3prime 70
    --samplesLabel some_sample another_one athird this_one_failed
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_scale():
    exp_npz = ROOT + "srmat_scale.gz"
    exp_mat = ROOT + "srmat_scale.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_IN1} {REGIONS_IN2} {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 98 -a 422 --binSize 2 -m 842
    --unscaled5prime 22 --unscaled3prime 44
    --scale 0.78
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_bed12():
    exp_npz = ROOT + "srmat_bed12.gz"
    exp_mat = ROOT + "srmat_bed12.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_BED12}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 98 -a 422 --binSize 2 -m 842
    --unscaled5prime 22 --unscaled3prime 44
    --scale 0.78
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_bed12mg():
    exp_npz = ROOT + "srmat_bed12mg.gz"
    exp_mat = ROOT + "srmat_bed12mg.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_BED12}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 98 -a 422 --binSize 2 -m 842
    --unscaled5prime 22 --unscaled3prime 44
    --scale 0.78
    --metagene
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"

def test_compute_matrix_scaleregions_gtfmg():
    exp_npz = ROOT + "srmat_gtfmg.gz"
    exp_mat = ROOT + "srmat_gtfmg.tab"

    _, outfile_npz = tempfile.mkstemp(suffix='.gz')
    _, outfile_mat = tempfile.mkstemp(suffix='.tab')
    args = f"""
    scale-regions
    --regionsFileName {REGIONS_GTF}
    --scoreFileName {BIGWIG_IN1} {BIGWIG_IN2} {BIGWIG_IN3} {BIGWIG_IN4}
    -b 98 -a 422 --binSize 2 -m 842
    --unscaled5prime 22 --unscaled3prime 44
    --scale 0.78
    --metagene
    -o {outfile_npz} --outFileNameMatrix {outfile_mat}
    """.split()
    cm.main(args)

    header_differences, data_differences, rowdiffdic = _compare_mat_gz(outfile_npz, exp_npz)
    assert not header_differences, f"header mismatch: {header_differences}"
    assert not data_differences, f"data mismatch: {data_differences}\nrowdiffdict: {rowdiffdic}"

    mat_match, mat_diffs = _compare_tab_files(outfile_mat, exp_mat)
    assert mat_match, f"matrix mismatch: {mat_diffs}"
