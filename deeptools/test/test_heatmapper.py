import os
import sys

import deeptools.computeMatrix
import deeptools.plotHeatmap
import deeptools.plotProfile
import deeptools.utilities
import json

__author__ = 'Fidel'

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_heatmapper/"


def cmpMatrices(f1, f2):
    """
    The header produced by computeMatrix will be different every time a command is run in python3!
    """
    rv = True
    file1 = open(f1)
    file2 = open(f2)
    for l1, l2 in zip(file1, file2):
        if isinstance(l1, bytes):
            l1 = l1.decode()
            l2 = l2.decode()
        l1 = l1.strip()
        l2 = l2.strip()
        if l1.startswith("@"):
            p1 = json.loads(l1[1:])
            p2 = json.loads(l2[1:])
            for k, v in p1.items():
                if k not in p2.keys():
                    sys.stderr.write("key in {} missing: {} not in {}\n".format(f1, k, p2.keys()))
                    rv = False
                if p1[k] != p2[k]:
                    sys.stderr.write("values of '{}' is different: {} not in {}\n".format(k, p1[k], p2[k]))
                    rv = False
            for k in p2.keys():
                if k not in p1.keys():
                    sys.stderr.write("key in {} missing: {} not in {}\n".format(f2, k, p1.keys()))
                    rv = False
        else:
            if l1 != l2:
                sys.stderr.write("lines differ:\n{}\n    vs\n{}\n".format(l1, l2))
                rv = False
    file1.close()
    file2.close()
    return rv


def test_computeMatrix_reference_point():
    args = "reference-point -R {0}/test2.bed -S {0}/test.bw  -b 100 -a 100 " \
           "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1".format(ROOT).split()
    deeptools.computeMatrix.main(args)
    os.system('gunzip -f /tmp/_test.mat.gz')
    assert cmpMatrices(ROOT + '/master.mat', '/tmp/_test.mat') is True
    os.remove('/tmp/_test.mat')


def test_computeMatrix_reference_point_center():
    args = "reference-point -R {0}/test2.bed -S {0}/test.bw  -b 100 -a 100 --referencePoint center " \
           "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1".format(ROOT).split()
    deeptools.computeMatrix.main(args)
    os.system('gunzip -f /tmp/_test.mat.gz')
    assert cmpMatrices(ROOT + '/master_center.mat', '/tmp/_test.mat') is True
    os.remove('/tmp/_test.mat')


def test_computeMatrix_reference_point_tes():
    args = "reference-point -R {0}/test2.bed -S {0}/test.bw  -b 100 -a 100 --referencePoint TES " \
           "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1".format(ROOT).split()
    deeptools.computeMatrix.main(args)
    os.system('gunzip -f /tmp/_test.mat.gz')
    assert cmpMatrices(ROOT + '/master_TES.mat', '/tmp/_test.mat') is True
    os.remove('/tmp/_test.mat')


def test_computeMatrix_reference_point_missing_data_as_zero():
    args = "reference-point -R {0}/test2.bed -S {0}/test.bw  -b 100 -a 100 " \
           "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1 --missingDataAsZero".format(ROOT).split()
    deeptools.computeMatrix.main(args)
    os.system('gunzip -f /tmp/_test.mat.gz')
    assert cmpMatrices(ROOT + '/master_nan_to_zero.mat', '/tmp/_test.mat') is True
    os.remove('/tmp/_test.mat')


def test_computeMatrix_scale_regions():
    args = "scale-regions -R {0}/test2.bed -S {0}/test.bw  -b 100 -a 100 -m 100 " \
           "--outFileName /tmp/_test2.mat.gz -bs 1 -p 1".format(ROOT).split()

    deeptools.computeMatrix.main(args)
    os.system('gunzip -f /tmp/_test2.mat.gz')
    assert cmpMatrices(ROOT + '/master_scale_reg.mat', '/tmp/_test2.mat') is True
    os.remove('/tmp/_test2.mat')


def test_computeMatrix_multiple_bed():
    args = "reference-point -R {0}/group1.bed {0}/group2.bed -S {0}/test.bw  -b 100 -a 100 " \
           "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1".format(ROOT).split()
    deeptools.computeMatrix.main(args)
    os.system('gunzip -f /tmp/_test.mat.gz')
    assert cmpMatrices(ROOT + '/master_multibed.mat', '/tmp/_test.mat') is True
    os.remove('/tmp/_test.mat')


def test_computeMatrix_region_extend_over_chr_end():
    args = "reference-point -R {0}/group1.bed {0}/group2.bed -S {0}/test.bw  -b 100 -a 500 " \
           "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1".format(ROOT).split()
    deeptools.computeMatrix.main(args)
    os.system('gunzip -f /tmp/_test.mat.gz')
    assert cmpMatrices(ROOT + '/master_extend_beyond_chr_size.mat', '/tmp/_test.mat') is True
    os.remove('/tmp/_test.mat')


def test_computeMatrix_unscaled():
    args = "scale-regions -S {0}/unscaled.bigWig -R {0}/unscaled.bed -a 300 -b 500 --unscaled5prime 100 --unscaled3prime 50 " \
           "--outFileName /tmp/_test.mat.gz -bs 10 -p 1".format(ROOT).split()
    deeptools.computeMatrix.main(args)
    os.system('gunzip -f /tmp/_test.mat.gz')
    assert cmpMatrices(ROOT + '/master_unscaled.mat', '/tmp/_test.mat') is True
    os.remove('/tmp/_test.mat')


def test_computeMatrix_gtf():
    args = "scale-regions -S {0}../test_data/test1.bw.bw -R {0}../test_data/test.gtf -a 300 -b 500 --unscaled5prime 20 --unscaled3prime 50 " \
           "--outFileName /tmp/_test_gtf.mat.gz -bs 10 -p 1".format(ROOT).split()
    deeptools.computeMatrix.main(args)
    os.system('gunzip -f /tmp/_test_gtf.mat.gz')
    assert cmpMatrices(ROOT + '/master_gtf.mat', '/tmp/_test_gtf.mat') is True
    os.remove('/tmp/_test_gtf.mat')


def test_computeMatrix_metagene():
    args = "scale-regions -S {0}../test_data/test1.bw.bw -R {0}../test_data/test.gtf -a 300 -b 500 --unscaled5prime 20 --unscaled3prime 50 " \
           "--outFileName /tmp/_test_metagene.mat.gz -bs 10 -p 1 --metagene".format(ROOT).split()
    deeptools.computeMatrix.main(args)
    os.system('gunzip -f /tmp/_test_metagene.mat.gz')
    assert cmpMatrices(ROOT + '/master_metagene.mat', '/tmp/_test_metagene.mat') is True
    os.remove('/tmp/_test_metagene.mat')


def test_chopRegions_body():
    region = [(0, 200), (300, 400), (800, 900)]
    lbins, bodybins, rbins, padLeft, padRight = deeptools.heatmapper.chopRegions(region, left=0, right=0)
    e_lbins = []
    e_rbins = []
    e_padLeft = 0
    e_padRight = 0
    assert f"{lbins}" == f"{e_lbins}"
    assert f"{rbins}" == f"{e_rbins}"
    assert f"{bodybins}" == f"{region}"
    assert f"{padLeft}" == f"{e_padLeft}"
    assert f"{padRight}" == f"{e_padRight}"
    # Unscaled 5', 3'
    lbins, bodybins, rbins, padLeft, padRight = deeptools.heatmapper.chopRegions(region, left=150, right=150)
    e_lbins = [(0, 150)]
    e_rbins = [(350, 400), (800, 900)]
    e_bodybins = [(150, 200), (300, 350)]
    e_padLeft = 0
    e_padRight = 0
    assert f"{lbins}" == f"{e_lbins}"
    assert f"{rbins}" == f"{e_rbins}"
    assert f"{bodybins}" == f"{e_bodybins}"
    assert f"{padLeft}" == f"{e_padLeft}"
    assert f"{padRight}" == f"{e_padRight}"


def test_chopRegions_TSS():
    region = [(0, 200), (300, 400), (800, 900)]
    # + strand, 250 downstream
    downstream, body, unscaled3prime, padRight, _ = deeptools.heatmapper.chopRegions(region, left=250)
    e_downstream = [(0, 200), (300, 350)]
    e_body = [(350, 400), (800, 900)]
    e_unscaled3prime = []
    e_padRight = 0
    e_ = 0
    assert f"{downstream}" == f"{e_downstream}"
    assert f"{body}" == f"{e_body}"
    assert f"{unscaled3prime}" == f"{e_unscaled3prime}"
    assert f"{padRight}" == f"{e_padRight}"
    assert f"{_}" == f"{e_}"
    # + strand, 500 downstream
    downstream, body, unscaled3prime, padRight, _ = deeptools.heatmapper.chopRegions(region, left=500)
    e_body = []
    e_unscaled3prime = []
    e_padRight = 100
    e_ = 0
    assert f"{downstream}" == f"{region}"
    assert f"{body}" == f"{e_body}"
    assert f"{unscaled3prime}" == f"{e_unscaled3prime}"
    assert f"{padRight}" == f"{e_padRight}"
    assert f"{_}" == f"{e_}"
    # - strand, 250 downstream (labeled "upstream" due to being on the - strand)
    unscaled5prime, body, upstream, _, padLeft = deeptools.heatmapper.chopRegions(region, right=250)
    e_upstream = [(150, 200), (300, 400), (800, 900)]
    e_body = [(0, 150)]
    e_unscaled5prime = []
    e_padLeft = 0
    e_ = 0
    assert f"{upstream}" == f"{e_upstream}"
    assert f"{body}" == f"{e_body}"
    assert f"{unscaled5prime}" == f"{e_unscaled5prime}"
    assert f"{padLeft}" == f"{e_padLeft}"
    assert f"{_}" == f"{e_}"
    # - strand, 500 downstream (labeled "upstream" due to being on the - strand)
    unscaled5prime, body, upstream, _, padLeft = deeptools.heatmapper.chopRegions(region, right=500)
    e_body = []
    e_unscaled5prime = []
    e_padLeft = 100
    e_ = 0
    assert f"{upstream}" == f"{region}"
    assert f"{body}" == f"{e_body}"
    assert f"{unscaled5prime}" == f"{e_unscaled5prime}"
    assert f"{padLeft}" == f"{e_padLeft}"
    assert f"{_}" == f"{e_}"


def test_chopRegions_TES():
    region = [(0, 200), (300, 400), (800, 900)]
    # + strand, 250 upstream
    unscaled5prime, body, upstream, _, padLeft = deeptools.heatmapper.chopRegions(region, right=250)
    e_unscaled5prime = []
    e_body = [(0, 150)]
    e_upstream = [(150, 200), (300, 400), (800, 900)]
    e_ = 0
    e_padLeft = 0
    assert f"{unscaled5prime}" == f"{e_unscaled5prime}"
    assert f"{body}" == f"{e_body}"
    assert f"{upstream}" == f"{e_upstream}"
    assert f"{_}" == f"{e_}"
    assert f"{padLeft}" == f"{e_padLeft}"
    # + strand, 500 upstream
    unscaled5prime, body, upstream, _, padLeft = deeptools.heatmapper.chopRegions(region, right=500)
    e_unscaled5prime = []
    e_body = []
    e_ = 0
    e_padLeft = 100
    assert f"{unscaled5prime}" == f"{e_unscaled5prime}"
    assert f"{body}" == f"{e_body}"
    assert f"{upstream}" == f"{region}"
    assert f"{_}" == f"{e_}"
    assert f"{padLeft}" == f"{e_padLeft}"
    # + strand, 250 downstream (labeled "upstream" due to being on the - strand)
    downstream, body, unscaled3prime, padRight, _ = deeptools.heatmapper.chopRegions(region, left=250)
    e_downstream = [(0, 200), (300, 350)]
    e_body = [(350, 400), (800, 900)]
    e_unscaled3prime = []
    e_padRight = 0
    e_ = 0
    assert f"{downstream}" == f"{e_downstream}"
    assert f"{body}" == f"{e_body}"
    assert f"{unscaled3prime}" == f"{e_unscaled3prime}"
    assert f"{padRight}" == f"{e_padRight}"
    assert f"{_}" == f"{e_}"
    # + strand, 500 downstream (labeled "upstream" due to being on the - strand)
    downstream, body, unscaled3prime, padRight, _ = deeptools.heatmapper.chopRegions(region, left=500)
    e_body = []
    e_unscaled3prime = []
    e_padRight = 100
    e_ = 0
    assert f"{downstream}" == f"{region}"
    assert f"{body}" == f"{e_body}"
    assert f"{unscaled3prime}" == f"{e_unscaled3prime}"
    assert f"{padRight}" == f"{e_padRight}"
    assert f"{_}" == f"{e_}"


def test_chopRegionsFromMiddle():
    region = [(0, 200), (300, 400), (800, 900)]
    # + strand, 100 upstream/200 downstream
    upstream, downstream, padLeft, padRight = deeptools.heatmapper.chopRegionsFromMiddle(region, left=100, right=200)
    e_upstream = [(100, 200)]
    e_downstream = [(300, 400), (800, 900)]
    e_padLeft = 0
    e_padRight = 0
    assert f"{upstream}" == f"{e_upstream}"
    assert f"{downstream}" == f"{e_downstream}"
    assert f"{padLeft}" == f"{e_padLeft}"
    assert f"{padRight}" == f"{e_padRight}"
    # + strand, 250 upstream/300 downstream
    upstream, downstream, padLeft, padRight = deeptools.heatmapper.chopRegionsFromMiddle(region, left=250, right=300)
    e_upstream = [(0, 200)]
    e_downstream = [(300, 400), (800, 900)]
    e_padLeft = 50
    e_padRight = 100
    assert f"{upstream}" == f"{e_upstream}"
    assert f"{downstream}" == f"{e_downstream}"
    assert f"{padLeft}" == f"{e_padLeft}"
    assert f"{padRight}" == f"{e_padRight}"
    # - strand, 100 upstream/200 downstream
    upstream, downstream, padLeft, padRight = deeptools.heatmapper.chopRegionsFromMiddle(region, left=200, right=100)
    e_upstream = [(0, 200)]
    e_downstream = [(300, 400)]
    e_padLeft = 0
    e_padRight = 0
    assert f"{upstream}" == f"{e_upstream}"
    assert f"{downstream}" == f"{e_downstream}"
    assert f"{padLeft}" == f"{e_padLeft}"
    assert f"{padRight}" == f"{e_padRight}"
    # - strand, 250 upstream/300 downstream
    upstream, downstream, padLeft, padRight = deeptools.heatmapper.chopRegionsFromMiddle(region, left=300, right=250)
    e_upstream = [(0, 200)]
    e_downstream = [(300, 400), (800, 900)]
    e_padLeft = 100
    e_padRight = 50
    assert f"{upstream}" == f"{e_upstream}"
    assert f"{downstream}" == f"{e_downstream}"
    assert f"{padLeft}" == f"{e_padLeft}"
    assert f"{padRight}" == f"{e_padRight}"
