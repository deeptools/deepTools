import deeptools.multiBamSummary2 as mbs
import numpy as np
import numpy.testing as nt
import math
import csv

import os.path
import tempfile

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/test_mbs/"
BAM = ROOT + "test1"
GTF = ROOT + "test.gtf"
BED = ROOT + "test.bed"
BED2 = ROOT + "test2.bed"
BED12 = ROOT + "bed12.bed"

BAMA = ROOT + "testA"
BAMB = ROOT + "testB"
BLACKLIST = ROOT + "blacklist.bed"
GTF2 = ROOT + "test_customids.gtf"

def compare_tsv(exp_tsv, obs_tsv, delta=1.0):
    with open(exp_tsv) as e, open(obs_tsv) as o:
        exp = list(csv.reader(e, delimiter="\t"))
        obs = list(csv.reader(o, delimiter="\t"))

    assert len(exp) == len(obs), f"row count differs: {len(exp)} != {len(obs)}"

    for i, (erow, orow) in enumerate(zip(exp, obs)):
        assert len(erow) == len(orow), f"col count differs on row {i}"
        for j, (a, b) in enumerate(zip(erow, orow)):
            try:
                assert math.isclose(float(a), float(b), abs_tol=delta), \
                    f"row {i}, col {j}: {a} != {b}"
            except ValueError:
                a = a.replace('.cram', '').replace('.bam', '')
                b = b.replace('.cram', '').replace('.bam', '')
                assert a == b, f"row {i}, col {j}: {a!r} != {b!r}"

def test_multiBamSummary_bedmode_gtf():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"BED-file --BED {GTF} -b {fname} {fname} -o {outfile}".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix, np.array([[144.0, 144.0],[143.0, 143.0]]))

def test_multiBamSummary_bedmode_bed():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"BED-file --BED {BED} -b {fname} {fname} -o {outfile}".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix,
            np.array(
                [[1.0, 1.0], [144.0, 144.0], [144.0, 144.0], [6.0, 6.0], [143.0, 143.0], [22.0, 22.0], [25.0, 25.0], [1.0, 1.0], [0.0, 0.0]]
            )
        )

def test_multiBamSummary_bedmode_multibed():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"BED-file --BED {BED} {GTF} {BED} -b {fname} {fname} -o {outfile}".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix,
            np.array(
                [[144.0, 144.0], [1.0, 1.0], [1.0, 1.0], [144.0, 144.0], [144.0, 144.0], [144.0, 144.0], [144.0, 144.0], [143.0, 143.0], [6.0, 6.0], [6.0, 6.0], [143.0, 143.0], [143.0, 143.0], [22.0, 22.0], [22.0, 22.0], [25.0, 25.0], [25.0, 25.0], [1.0, 1.0], [1.0, 1.0], [0.0, 0.0], [0.0, 0.0]]
            )
        )

def test_multiBamSummary_bedmode_diffbam():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname1 = BAMA + fname
        fname2 = BAMB + fname
        args = f"BED-file --BED {BED} -b {fname1} {fname2} -o {outfile}".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix,
            np.array(
                [[0.0, 0.0], [2.0, 4.0], [2.0, 4.0], [1.0, 1.0], [2.0, 3.0]]
            )
        )

def test_multiBamSummary_bedmode_labels():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname1 = BAMA + fname
        fname2 = BAMB + fname
        args = f"BED-file --BED {BED} -b {fname1} {fname2} -o {outfile} --labels bam_part1 bam2".split()
        mbs.main(args)
        resp = np.load(outfile)
        labels = resp['labels']
        decoded_labels = [bytes(lab).decode('ascii').rstrip('\x00') for lab in labels]

        assert decoded_labels == ['bam_part1', 'bam2']

def test_multiBamSummary_bedmode_smartlabels():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname1 = BAMA + fname
        fname2 = BAMB + fname
        args = f"BED-file --BED {BED} -b {fname1} {fname2} -o {outfile} --smartLabels".split()
        mbs.main(args)
        resp = np.load(outfile)
        labels = resp['labels']
        decoded_labels = [bytes(lab).decode('ascii').rstrip('\x00') for lab in labels]

        assert decoded_labels == ['testA', 'testB']

def test_multiBamSummary_bedmode_region():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"BED-file --BED {BED} -b {fname} {fname} -o {outfile} --region 3R:500:1000".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix,
            np.array(
                [[144.0, 144.0], [144.0, 144.0], [143.0, 143.0], [22.0, 22.0], [25.0, 25.0], [1.0, 1.0]]
            )
        )

def test_multiBamSummary_bedmode_blacklist():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"BED-file --BED {BED} -b {fname} {fname} -o {outfile} --blackListFileName {BLACKLIST}".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix,
            np.array(
                [[1.0, 1.0], [59.0, 59.0], [59.0, 59.0], [6.0, 6.0], [58.0, 58.0], [16.0, 16.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
            )
        )

def test_multiBamSummary_bedmode_outcounts_sf():
    exp_count = ROOT + "bedmode_counts.tsv"
    exp_sf = ROOT + "bedmode_sf.txt"

    _, outfile = tempfile.mkstemp(suffix=".npz")
    _, osf = tempfile.mkstemp(suffix=".txt")
    _, ocount = tempfile.mkstemp(suffix=".tsv")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"BED-file --BED {BED} -b {fname} {fname} -o {outfile} --outRawCounts {ocount} --scalingFactors {osf}".split()
        mbs.main(args)
        compare_tsv(exp_count, ocount)
        compare_tsv(exp_sf, osf)

def test_multiBamSummary_bedmode_outcounts_sf2():
    exp_count = ROOT + "bedmode_counts2.tsv"
    exp_sf = ROOT + "bedmode_sf2.txt"

    _, outfile = tempfile.mkstemp(suffix=".npz")
    _, osf = tempfile.mkstemp(suffix=".txt")
    _, ocount = tempfile.mkstemp(suffix=".tsv")
    for fname in ['.bam', '.cram']:
        fname1 = BAMA + fname
        fname2 = BAMB + fname
        args = f"BED-file --BED {BED} -b {fname1} {fname2} -o {outfile} --outRawCounts {ocount} --scalingFactors {osf}".split()
        mbs.main(args)
        compare_tsv(exp_count, ocount)
        compare_tsv(exp_sf, osf)

def test_multiBamSummary_bedmode_center_ext():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"BED-file --BED {BED} -b {fname} {fname} -o {outfile} -e 40 --centerReads".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix,
            np.array(
                [[1.0, 1.0], [144.0, 144.0], [144.0, 144.0], [6.0, 6.0], [143.0, 143.0], [22.0, 22.0], [25.0, 25.0], [1.0, 1.0], [0.0, 0.0]]
            )
        )

def test_multiBamSummary_bedmode_center_ext2():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"BED-file --BED {BED} -b {fname} {fname} -o {outfile} -e 80 --centerReads".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix,
            np.array(
                [[1.0, 1.0], [144.0, 144.0], [144.0, 144.0], [7.0, 7.0], [143.0, 143.0], [23.0, 23.0], [22.0, 22.0], [3.0, 3.0], [0.0, 0.0]]
            )
        )

def test_multiBamSummary_bedmode_center_ext3():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"BED-file --BED {BED} -b {fname} {fname} -o {outfile} -e 150".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix,
            np.array(
                [[3.0, 3.0], [144.0, 144.0], [144.0, 144.0], [16.0, 16.0], [144.0, 144.0], [35.0, 35.0], [41.0, 41.0], [16.0, 16.0], [0.0, 0.0]]
            )
        )

def test_multiBamSummary_bedmode_center_ext4():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"BED-file --BED {BED} -b {fname} {fname} -o {outfile} --centerReads -e 150".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix,
            np.array(
                [[1.0, 1.0], [144.0, 144.0], [144.0, 144.0], [10.0, 10.0], [143.0, 143.0], [23.0, 23.0], [23.0, 23.0], [5.0, 5.0], [0.0, 0.0]]
            )
        )

def test_multiBamSummary_gtf_customids():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = (
            f"BED-file --BED {GTF2} -b {fname} {fname} -o {outfile} "
            "--transcriptID mRNA --exonID CDS --transcript_id_designator ID"
        ).split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix, np.array([[144.0, 144.0], [143.0, 143.0]]))

def test_multiBamSummary_gtf_customids_metagene():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = (
            f"BED-file --BED {GTF2} -b {fname} {fname} -o {outfile} --metagene "
            "--transcriptID mRNA --exonID CDS --transcript_id_designator ID"
        ).split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']

        nt.assert_allclose(matrix, np.array([[25.0, 25.0], [31.0, 31.0]]))

def test_multiBamSummary_bedmode_readfilters():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = (
            f"BED-file --BED {BED2} -b {fname} {fname} -o {outfile} "
            "--samFlagInclude 16 --samFlagExclude 1024 "
            "--minFragmentLength 40 --maxFragmentLength 60 "
            "--minMappingQuality 20"
        ).split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']
        nt.assert_allclose(matrix, np.array([[36.0, 36.0], [36.0, 36.0]]))

def test_multiBamSummary_bedmode_gtfmg():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = (
            f"BED-file --BED {GTF} -b {fname} {fname} -o {outfile} "
            "--metagene "
        ).split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']
        nt.assert_allclose(matrix, np.array([[25.0, 25.0], [31.0, 31.0]]))

def test_multiBamSummary_bedmode_bed12mg():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = (
            f"BED-file --BED {BED12} -b {fname} {fname} -o {outfile} "
            "--metagene "
        ).split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']
        nt.assert_allclose(matrix, np.array([[25.0, 25.0], [31.0, 31.0]]))

def test_multiBamSummary_binmode():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = (
            f"bins --binSize 10 -b {fname} {fname} -o {outfile}"
        ).split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']
        nt.assert_allclose(matrix, np.array(
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [3.0, 3.0], [4.0, 4.0], [4.0, 4.0], [4.0, 4.0], [4.0, 4.0], [7.0, 7.0], [10.0, 10.0], [17.0, 17.0], [21.0, 21.0], [21.0, 21.0], [20.0, 20.0], [17.0, 17.0], [17.0, 17.0], [12.0, 12.0], [8.0, 8.0], [7.0, 7.0], [8.0, 8.0], [8.0, 8.0], [7.0, 7.0], [5.0, 5.0], [6.0, 6.0], [6.0, 6.0], [5.0, 5.0], [4.0, 4.0], [5.0, 5.0], [4.0, 4.0], [3.0, 3.0], [4.0, 4.0], [6.0, 6.0], [8.0, 8.0], [8.0, 8.0], [10.0, 10.0], [12.0, 12.0], [13.0, 13.0], [13.0, 13.0], [12.0, 12.0], [11.0, 11.0], [10.0, 10.0], [9.0, 9.0], [7.0, 7.0], [6.0, 6.0], [5.0, 5.0], [5.0, 5.0], [5.0, 5.0], [6.0, 6.0], [6.0, 6.0], [5.0, 5.0], [9.0, 9.0], [10.0, 10.0], [10.0, 10.0], [14.0, 14.0], [14.0, 14.0], [14.0, 14.0], [11.0, 11.0], [8.0, 8.0], [10.0, 10.0], [5.0, 5.0], [6.0, 6.0], [8.0, 8.0], [8.0, 8.0], [11.0, 11.0], [9.0, 9.0], [15.0, 15.0], [14.0, 14.0], [12.0, 12.0], [13.0, 13.0], [12.0, 12.0], [13.0, 13.0], [6.0, 6.0], [6.0, 6.0], [6.0, 6.0], [4.0, 4.0], [3.0, 3.0], [1.0, 1.0], [3.0, 3.0], [5.0, 5.0], [15.0, 15.0], [22.0, 22.0], [23.0, 23.0], [26.0, 26.0], [25.0, 25.0], [26.0, 26.0], [18.0, 18.0], [11.0, 11.0], [9.0, 9.0], [6.0, 6.0], [6.0, 6.0], [3.0, 3.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        ))

def test_multiBamSummary_binmode2():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = (
            f"bins --binSize 18 -b {fname} {fname} -o {outfile}"
        ).split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']
        nt.assert_allclose(matrix, np.array(
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [4.0, 4.0], [4.0, 4.0], [8.0, 8.0], [18.0, 18.0], [22.0, 22.0], [21.0, 21.0], [19.0, 19.0], [10.0, 10.0], [8.0, 8.0], [9.0, 9.0], [6.0, 6.0], [6.0, 6.0], [8.0, 8.0], [5.0, 5.0], [4.0, 4.0], [8.0, 8.0], [9.0, 9.0], [12.0, 12.0], [15.0, 15.0], [14.0, 14.0], [11.0, 11.0], [7.0, 7.0], [6.0, 6.0], [5.0, 5.0], [6.0, 6.0], [8.0, 8.0], [12.0, 12.0], [14.0, 14.0], [14.0, 14.0], [11.0, 11.0], [11.0, 11.0], [6.0, 6.0], [9.0, 9.0], [12.0, 12.0], [15.0, 15.0], [15.0, 15.0], [14.0, 14.0], [12.0, 12.0], [6.0, 6.0], [5.0, 5.0], [3.0, 3.0], [5.0, 5.0], [22.0, 22.0], [26.0, 26.0], [27.0, 27.0], [18.0, 18.0], [9.0, 9.0], [6.0, 6.0], [3.0, 3.0], [1.0, 1.0], [1.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        ))

def test_multiBamSummary_binmode3():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = (
            f"bins --binSize 18 --distanceBetweenBins 3 -b {fname} {fname} -o {outfile}"
        ).split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']
        nt.assert_allclose(matrix, np.array(
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [3.0, 3.0], [4.0, 4.0], [4.0, 4.0], [8.0, 8.0], [20.0, 20.0], [21.0, 21.0], [19.0, 19.0], [12.0, 12.0], [8.0, 8.0], [9.0, 9.0], [6.0, 6.0], [6.0, 6.0], [6.0, 6.0], [3.0, 3.0], [8.0, 8.0], [9.0, 9.0], [14.0, 14.0], [14.0, 14.0], [12.0, 12.0], [9.0, 9.0], [6.0, 6.0], [6.0, 6.0], [6.0, 6.0], [10.0, 10.0], [14.0, 14.0], [14.0, 14.0], [11.0, 11.0], [11.0, 11.0], [8.0, 8.0], [12.0, 12.0], [15.0, 15.0], [14.0, 14.0], [14.0, 14.0], [8.0, 8.0], [6.0, 6.0], [3.0, 3.0], [3.0, 3.0], [21.0, 21.0], [26.0, 26.0], [28.0, 28.0], [17.0, 17.0], [7.0, 7.0], [4.0, 4.0], [1.0, 1.0], [1.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        ))

def test_multiBamSummary_binmode_labels():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname1 = BAMA + fname
        fname2 = BAMB + fname
        args = f"bins --binSize 10 -b {fname1} {fname2} -o {outfile} --labels bam_part1 bam2".split()
        mbs.main(args)
        resp = np.load(outfile)
        labels = resp['labels']
        decoded_labels = [bytes(lab).decode('ascii').rstrip('\x00') for lab in labels]

        assert decoded_labels == ['bam_part1', 'bam2']

def test_multiBamSummary_binmode_smartlabels():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname1 = BAMA + fname
        fname2 = BAMB + fname
        args = f"bins --binSize 10 -b {fname1} {fname2} -o {outfile} --smartLabels".split()
        mbs.main(args)
        resp = np.load(outfile)
        labels = resp['labels']
        decoded_labels = [bytes(lab).decode('ascii').rstrip('\x00') for lab in labels]

        assert decoded_labels == ['testA', 'testB']

def test_multiBamSummary_binmode_region():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"bins --binSize 10 -b {fname} {fname} -o {outfile} --region 3R:500:1000 ".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']
        nt.assert_allclose(matrix, np.array(
            [[1.0, 1.0], [1.0, 1.0], [3.0, 3.0], [4.0, 4.0], [6.0, 6.0], [6.0, 6.0], [5.0, 5.0], [9.0, 9.0], [10.0, 10.0], [10.0, 10.0], [14.0, 14.0], [14.0, 14.0], [14.0, 14.0], [11.0, 11.0], [8.0, 8.0], [10.0, 10.0], [5.0, 5.0], [6.0, 6.0], [8.0, 8.0], [8.0, 8.0], [11.0, 11.0], [9.0, 9.0], [15.0, 15.0], [14.0, 14.0], [12.0, 12.0], [13.0, 13.0], [12.0, 12.0], [13.0, 13.0], [6.0, 6.0], [6.0, 6.0], [6.0, 6.0], [4.0, 4.0], [3.0, 3.0], [1.0, 1.0], [3.0, 3.0], [5.0, 5.0], [15.0, 15.0], [22.0, 22.0], [23.0, 23.0], [26.0, 26.0], [25.0, 25.0], [26.0, 26.0], [18.0, 18.0], [11.0, 11.0], [9.0, 9.0], [6.0, 6.0], [5.0, 5.0], [2.0, 2.0], [0.0, 0.0], [0.0, 0.0]]
        ))

def test_multiBamSummary_binmode_blacklist():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"bins --binSize 10 -b {fname} {fname} -o {outfile} --blackListFileName {BLACKLIST} ".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']
        nt.assert_allclose(matrix, np.array(
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [3.0, 3.0], [4.0, 4.0], [4.0, 4.0], [4.0, 4.0], [4.0, 4.0], [7.0, 7.0], [10.0, 10.0], [17.0, 17.0], [21.0, 21.0], [21.0, 21.0], [20.0, 20.0], [17.0, 17.0], [17.0, 17.0], [12.0, 12.0], [8.0, 8.0], [7.0, 7.0], [8.0, 8.0], [8.0, 8.0], [7.0, 7.0], [5.0, 5.0], [6.0, 6.0], [6.0, 6.0], [5.0, 5.0], [4.0, 4.0], [5.0, 5.0], [4.0, 4.0], [3.0, 3.0], [4.0, 4.0], [6.0, 6.0], [8.0, 8.0], [8.0, 8.0], [10.0, 10.0], [12.0, 12.0], [13.0, 13.0], [13.0, 13.0], [11.0, 11.0], [8.0, 8.0], [6.0, 6.0], [4.0, 4.0], [2.0, 2.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        ))

def test_multiBamSummary_binmode_outcounts_sf():
    exp_count = ROOT + "binmode_counts.tsv"
    exp_sf = ROOT + "binmode_sf.txt"

    _, outfile = tempfile.mkstemp(suffix=".npz")
    _, osf = tempfile.mkstemp(suffix=".txt")
    _, ocount = tempfile.mkstemp(suffix=".tsv")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"bins --binSize 10 -b {fname} {fname} -o {outfile} --outRawCounts {ocount} --scalingFactors {osf}".split()
        mbs.main(args)
        compare_tsv(exp_count, ocount)
        compare_tsv(exp_sf, osf)

def test_multiBamSummary_binmode_outcounts_sf2():
    exp_count = ROOT + "binmode_counts2.tsv"
    exp_sf = ROOT + "binmode_sf2.txt"

    _, outfile = tempfile.mkstemp(suffix=".npz")
    _, osf = tempfile.mkstemp(suffix=".txt")
    _, ocount = tempfile.mkstemp(suffix=".tsv")
    for fname in ['.bam', '.cram']:
        fname1 = BAMA + fname
        fname2 = BAMB + fname
        args = f"bins --binSize 10 -b {fname1} {fname2} -o {outfile} --outRawCounts {ocount} --scalingFactors {osf}".split()
        mbs.main(args)
        compare_tsv(exp_count, ocount)
        compare_tsv(exp_sf, osf)

def test_multiBamSummary_binmode_center_ext():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"bins --binSize 10 -b {fname} {fname} -o {outfile} --centerReads --extendReads 20 ".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']
        nt.assert_allclose(matrix, np.array(
            [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [3.0, 3.0], [4.0, 4.0], [4.0, 4.0], [4.0, 4.0], [4.0, 4.0], [7.0, 7.0], [10.0, 10.0], [17.0, 17.0], [21.0, 21.0], [21.0, 21.0], [20.0, 20.0], [17.0, 17.0], [17.0, 17.0], [12.0, 12.0], [8.0, 8.0], [7.0, 7.0], [8.0, 8.0], [8.0, 8.0], [7.0, 7.0], [5.0, 5.0], [6.0, 6.0], [6.0, 6.0], [5.0, 5.0], [4.0, 4.0], [5.0, 5.0], [4.0, 4.0], [3.0, 3.0], [4.0, 4.0], [6.0, 6.0], [8.0, 8.0], [8.0, 8.0], [10.0, 10.0], [12.0, 12.0], [13.0, 13.0], [13.0, 13.0], [12.0, 12.0], [11.0, 11.0], [10.0, 10.0], [9.0, 9.0], [7.0, 7.0], [6.0, 6.0], [5.0, 5.0], [5.0, 5.0], [5.0, 5.0], [6.0, 6.0], [6.0, 6.0], [5.0, 5.0], [9.0, 9.0], [10.0, 10.0], [10.0, 10.0], [14.0, 14.0], [14.0, 14.0], [14.0, 14.0], [11.0, 11.0], [8.0, 8.0], [10.0, 10.0], [5.0, 5.0], [6.0, 6.0], [8.0, 8.0], [8.0, 8.0], [11.0, 11.0], [9.0, 9.0], [15.0, 15.0], [14.0, 14.0], [12.0, 12.0], [13.0, 13.0], [12.0, 12.0], [13.0, 13.0], [6.0, 6.0], [6.0, 6.0], [6.0, 6.0], [4.0, 4.0], [3.0, 3.0], [1.0, 1.0], [3.0, 3.0], [5.0, 5.0], [15.0, 15.0], [22.0, 22.0], [23.0, 23.0], [26.0, 26.0], [25.0, 25.0], [26.0, 26.0], [18.0, 18.0], [11.0, 11.0], [9.0, 9.0], [6.0, 6.0], [6.0, 6.0], [3.0, 3.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        ))

def test_multiBamSummary_binmode_center_ext2():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"bins --binSize 10 -b {fname} {fname} -o {outfile} --centerReads --extendReads 80 ".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']
        nt.assert_allclose(matrix, np.array(
            [[0.0, 0.0], [0.0, 0.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [4.0, 4.0], [6.0, 6.0], [7.0, 7.0], [10.0, 10.0], [11.0, 11.0], [18.0, 18.0], [20.0, 20.0], [19.0, 19.0], [18.0, 18.0], [16.0, 16.0], [16.0, 16.0], [11.0, 11.0], [9.0, 9.0], [8.0, 8.0], [9.0, 9.0], [7.0, 7.0], [7.0, 7.0], [5.0, 5.0], [5.0, 5.0], [8.0, 8.0], [6.0, 6.0], [6.0, 6.0], [4.0, 4.0], [4.0, 4.0], [4.0, 4.0], [4.0, 4.0], [7.0, 7.0], [8.0, 8.0], [10.0, 10.0], [10.0, 10.0], [11.0, 11.0], [11.0, 11.0], [9.0, 9.0], [10.0, 10.0], [10.0, 10.0], [11.0, 11.0], [10.0, 10.0], [8.0, 8.0], [8.0, 8.0], [8.0, 8.0], [6.0, 6.0], [5.0, 5.0], [7.0, 7.0], [7.0, 7.0], [7.0, 7.0], [6.0, 6.0], [9.0, 9.0], [11.0, 11.0], [11.0, 11.0], [12.0, 12.0], [12.0, 12.0], [12.0, 12.0], [9.0, 9.0], [9.0, 9.0], [9.0, 9.0], [7.0, 7.0], [6.0, 6.0], [5.0, 5.0], [9.0, 9.0], [10.0, 10.0], [11.0, 11.0], [13.0, 13.0], [15.0, 15.0], [16.0, 16.0], [12.0, 12.0], [11.0, 11.0], [7.0, 7.0], [5.0, 5.0], [3.0, 3.0], [2.0, 2.0], [3.0, 3.0], [1.0, 1.0], [12.0, 12.0], [16.0, 16.0], [20.0, 20.0], [25.0, 25.0], [26.0, 26.0], [28.0, 28.0], [18.0, 18.0], [15.0, 15.0], [12.0, 12.0], [8.0, 8.0], [7.0, 7.0], [5.0, 5.0], [4.0, 4.0], [3.0, 3.0], [2.0, 2.0], [1.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        ))

def test_multiBamSummary_binmode_center_ext3():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"bins --binSize 10 -b {fname} {fname} -o {outfile} --extendReads 150 ".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']
        nt.assert_allclose(matrix, np.array(
            [[2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [5.0, 5.0], [7.0, 7.0], [9.0, 9.0], [11.0, 11.0], [11.0, 11.0], [12.0, 12.0], [13.0, 13.0], [14.0, 14.0], [15.0, 15.0], [18.0, 18.0], [19.0, 19.0], [26.0, 26.0], [30.0, 30.0], [30.0, 30.0], [30.0, 30.0], [30.0, 30.0], [31.0, 31.0], [30.0, 30.0], [31.0, 31.0], [29.0, 29.0], [29.0, 29.0], [29.0, 29.0], [29.0, 29.0], [29.0, 29.0], [31.0, 31.0], [29.0, 29.0], [25.0, 25.0], [21.0, 21.0], [17.0, 17.0], [17.0, 17.0], [17.0, 17.0], [18.0, 18.0], [18.0, 18.0], [18.0, 18.0], [17.0, 17.0], [18.0, 18.0], [21.0, 21.0], [23.0, 23.0], [24.0, 24.0], [23.0, 23.0], [22.0, 22.0], [24.0, 24.0], [23.0, 23.0], [24.0, 24.0], [26.0, 26.0], [25.0, 25.0], [25.0, 25.0], [24.0, 24.0], [23.0, 23.0], [23.0, 23.0], [22.0, 22.0], [23.0, 23.0], [25.0, 25.0], [22.0, 22.0], [23.0, 23.0], [23.0, 23.0], [23.0, 23.0], [20.0, 20.0], [20.0, 20.0], [21.0, 21.0], [22.0, 22.0], [22.0, 22.0], [22.0, 22.0], [22.0, 22.0], [24.0, 24.0], [24.0, 24.0], [28.0, 28.0], [26.0, 26.0], [21.0, 21.0], [22.0, 22.0], [27.0, 27.0], [31.0, 31.0], [31.0, 31.0], [32.0, 32.0], [33.0, 33.0], [32.0, 32.0], [31.0, 31.0], [31.0, 31.0], [33.0, 33.0], [34.0, 34.0], [31.0, 31.0], [37.0, 37.0], [32.0, 32.0], [34.0, 34.0], [34.0, 34.0], [34.0, 34.0], [27.0, 27.0], [23.0, 23.0], [21.0, 21.0], [20.0, 20.0], [19.0, 19.0], [17.0, 17.0], [16.0, 16.0], [15.0, 15.0], [14.0, 14.0], [12.0, 12.0], [10.0, 10.0], [5.0, 5.0], [4.0, 4.0], [3.0, 3.0], [3.0, 3.0], [2.0, 2.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        ))

def test_multiBamSummary_binmode_center_ext4():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = f"bins --binSize 10 -b {fname} {fname} -o {outfile} --centerReads --extendReads 150 ".split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']
        nt.assert_allclose(matrix, np.array(
            [[0.0, 0.0], [0.0, 0.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [3.0, 3.0], [3.0, 3.0], [5.0, 5.0], [7.0, 7.0], [8.0, 8.0], [8.0, 8.0], [8.0, 8.0], [8.0, 8.0], [8.0, 8.0], [7.0, 7.0], [8.0, 8.0], [13.0, 13.0], [16.0, 16.0], [19.0, 19.0], [18.0, 18.0], [17.0, 17.0], [15.0, 15.0], [12.0, 12.0], [9.0, 9.0], [8.0, 8.0], [7.0, 7.0], [7.0, 7.0], [7.0, 7.0], [7.0, 7.0], [7.0, 7.0], [7.0, 7.0], [8.0, 8.0], [9.0, 9.0], [8.0, 8.0], [6.0, 6.0], [6.0, 6.0], [3.0, 3.0], [4.0, 4.0], [5.0, 5.0], [6.0, 6.0], [8.0, 8.0], [8.0, 8.0], [11.0, 11.0], [12.0, 12.0], [13.0, 13.0], [12.0, 12.0], [12.0, 12.0], [14.0, 14.0], [11.0, 11.0], [9.0, 9.0], [8.0, 8.0], [8.0, 8.0], [6.0, 6.0], [4.0, 4.0], [5.0, 5.0], [5.0, 5.0], [4.0, 4.0], [6.0, 6.0], [11.0, 11.0], [10.0, 10.0], [13.0, 13.0], [12.0, 12.0], [13.0, 13.0], [11.0, 11.0], [6.0, 6.0], [9.0, 9.0], [7.0, 7.0], [7.0, 7.0], [5.0, 5.0], [6.0, 6.0], [10.0, 10.0], [7.0, 7.0], [10.0, 10.0], [10.0, 10.0], [10.0, 10.0], [10.0, 10.0], [15.0, 15.0], [19.0, 19.0], [16.0, 16.0], [17.0, 17.0], [18.0, 18.0], [19.0, 19.0], [11.0, 11.0], [8.0, 8.0], [7.0, 7.0], [8.0, 8.0], [9.0, 9.0], [12.0, 12.0], [12.0, 12.0], [12.0, 12.0], [11.0, 11.0], [10.0, 10.0], [10.0, 10.0], [5.0, 5.0], [4.0, 4.0], [3.0, 3.0], [3.0, 3.0], [2.0, 2.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        ))

def test_multiBamSummary_binmode_readfilters():
    _, outfile = tempfile.mkstemp(suffix=".npz")
    for fname in ['.bam', '.cram']:
        fname = BAM + fname
        args = (
            f"bins --binSize 5 -b {fname} {fname} -o {outfile} "
            "--samFlagInclude 16 --samFlagExclude 1024 "
            "--minFragmentLength 40 --maxFragmentLength 60 "
            "--minMappingQuality 20"
        ).split()
        mbs.main(args)
        resp = np.load(outfile)
        matrix = resp['matrix']
        nt.assert_allclose(matrix, np.array(
            [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [4.0, 4.0], [4.0, 4.0], [4.0, 4.0], [4.0, 4.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [5.0, 5.0], [4.0, 4.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [2.0, 2.0], [3.0, 3.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [3.0, 3.0], [3.0, 3.0], [2.0, 2.0], [4.0, 4.0], [4.0, 4.0], [5.0, 5.0], [5.0, 5.0], [5.0, 5.0], [5.0, 5.0], [5.0, 5.0], [5.0, 5.0], [4.0, 4.0], [4.0, 4.0], [3.0, 3.0], [1.0, 1.0], [1.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [2.0, 2.0], [2.0, 2.0], [3.0, 3.0], [3.0, 3.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [1.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [5.0, 5.0], [5.0, 5.0], [7.0, 7.0], [8.0, 8.0], [8.0, 8.0], [9.0, 9.0], [9.0, 9.0], [10.0, 10.0], [10.0, 10.0], [10.0, 10.0], [11.0, 11.0], [6.0, 6.0], [6.0, 6.0], [4.0, 4.0], [3.0, 3.0], [3.0, 3.0], [2.0, 2.0], [2.0, 2.0], [1.0, 1.0], [1.0, 1.0], [2.0, 2.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        ))
