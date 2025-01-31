import deeptools.computeMatrix as cm

import os.path
from os import unlink

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
REGIONS_IN1 = ROOT + "computeMatrix1.bed"
REGIONS_IN2 = ROOT + "computeMatrix2.bed"
BIGWIG_IN1 = ROOT + "bamCoverage_result4.bw"
BIGWIG_IN2 = ROOT + "computeMatrix2.bw"
OUT_ARCHIEVE1 = ROOT + "computeMatrix_result1.gz"
OUT_ARCHIEVE2 = ROOT + "computeMatrix_result2.gz"
OUT_ARCHIEVE3 = ROOT + "computeMatrix_result3.gz"

def test_compute_matrix_with_reference_point_and_advance_options_1():
    """
    Test minimal command line args compute matrix with reference based
    mode along with advance options of sorting using sum and average 
    type bin as sum
    """
    outfile = '/tmp/computematrix_1.gz'
    args = "reference-point --regionsFileName {} --scoreFileName {} -o {} -bs 10 --sortUsing sum --averageTypeBins sum -b 10 -a 10".format(REGIONS_IN1, BIGWIG_IN1, outfile).split()
    cm.main(args)

    archieve_file_size = os.path.getsize(OUT_ARCHIEVE1)
    expected_file_size = os.path.getsize(outfile)
    size_tolerance = 500
    size_difference = abs(archieve_file_size - expected_file_size)
    assert size_difference <= size_tolerance, "File size do not match"
    unlink(outfile)

def test_compute_matrix_with_reference_point_and_advance_options_2():
    """
    Test minimal command line args compute matrix with reference based mode 
    with before and after region start length
    """
    outfile = '/tmp/computematrix_2.gz'
    args = "reference-point --regionsFileName {} --scoreFileName {} -o {} -bs 10 -b 10 -a 10".format(REGIONS_IN2, BIGWIG_IN2, outfile).split()
    cm.main(args)

    archieve_file_size = os.path.getsize(OUT_ARCHIEVE2)
    expected_file_size = os.path.getsize(outfile)
    size_tolerance = 500
    size_difference = abs(archieve_file_size - expected_file_size)
    assert size_difference <= size_tolerance, "File size do not match"
    unlink(outfile)

def test_compute_matrix_with_scale_regions():
    """
    Test minimal command line args compute matrix with scale regions mode
    """
    outfile = '/tmp/computematrix_3.gz'
    args = "scale-regions --regionsFileName {} --scoreFileName {} -o {}".format(REGIONS_IN2, BIGWIG_IN2, outfile).split()
    cm.main(args)

    archieve_file_size = os.path.getsize(OUT_ARCHIEVE3)
    expected_file_size = os.path.getsize(outfile)
    size_tolerance = 500
    size_difference = abs(archieve_file_size - expected_file_size)
    assert size_difference <= size_tolerance, "File size do not match"
    unlink(outfile)