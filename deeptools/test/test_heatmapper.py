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
    assert(lbins == [])
    assert(rbins == [])
    assert(bodybins == region)
    assert(padLeft == 0)
    assert(padRight == 0)
    # Unscaled 5', 3'
    lbins, bodybins, rbins, padLeft, padRight = deeptools.heatmapper.chopRegions(region, left=150, right=150)
    assert(lbins == [(0, 150)])
    assert(rbins == [(350, 400), (800, 900)])
    assert(bodybins == [(150, 200), (300, 350)])
    assert(padLeft == 0)
    assert(padRight == 0)


def test_chopRegions_TSS():
    region = [(0, 200), (300, 400), (800, 900)]
    # + strand, 250 downstream
    downstream, body, unscaled3prime, padRight, _ = deeptools.heatmapper.chopRegions(region, left=250)
    assert(downstream == [(0, 200), (300, 350)])
    assert(body == [(350, 400), (800, 900)])
    assert(unscaled3prime == [])
    assert(padRight == 0)
    assert(_ == 0)
    # + strand, 500 downstream
    downstream, body, unscaled3prime, padRight, _ = deeptools.heatmapper.chopRegions(region, left=500)
    assert(downstream == region)
    assert(body == [])
    assert(unscaled3prime == [])
    assert(padRight == 100)
    assert(_ == 0)
    # - strand, 250 downstream (labeled "upstream" due to being on the - strand)
    unscaled5prime, body, upstream, _, padLeft = deeptools.heatmapper.chopRegions(region, right=250)
    assert(upstream == [(150, 200), (300, 400), (800, 900)])
    assert(body == [(0, 150)])
    assert(unscaled5prime == [])
    assert(padLeft == 0)
    assert(_ == 0)
    # - strand, 500 downstream (labeled "upstream" due to being on the - strand)
    unscaled5prime, body, upstream, _, padLeft = deeptools.heatmapper.chopRegions(region, right=500)
    assert(upstream == region)
    assert(body == [])
    assert(unscaled5prime == [])
    assert(padLeft == 100)
    assert(_ == 0)


def test_chopRegions_TES():
    region = [(0, 200), (300, 400), (800, 900)]
    # + strand, 250 upstream
    unscaled5prime, body, upstream, _, padLeft = deeptools.heatmapper.chopRegions(region, right=250)
    assert(unscaled5prime == [])
    assert(body == [(0, 150)])
    assert(upstream == [(150, 200), (300, 400), (800, 900)])
    assert(_ == 0)
    assert(padLeft == 0)
    # + strand, 500 upstream
    unscaled5prime, body, upstream, _, padLeft = deeptools.heatmapper.chopRegions(region, right=500)
    assert(unscaled5prime == [])
    assert(body == [])
    assert(upstream == region)
    assert(_ == 0)
    assert(padLeft == 100)
    # + strand, 250 downstream (labeled "upstream" due to being on the - strand)
    downstream, body, unscaled3prime, padRight, _ = deeptools.heatmapper.chopRegions(region, left=250)
    assert(downstream == [(0, 200), (300, 350)])
    assert(body == [(350, 400), (800, 900)])
    assert(unscaled3prime == [])
    assert(padRight == 0)
    assert(_ == 0)
    # + strand, 500 downstream (labeled "upstream" due to being on the - strand)
    downstream, body, unscaled3prime, padRight, _ = deeptools.heatmapper.chopRegions(region, left=500)
    assert(downstream == region)
    assert(body == [])
    assert(unscaled3prime == [])
    assert(padRight == 100)
    assert(_ == 0)


def test_chopRegionsFromMiddle():
    region = [(0, 200), (300, 400), (800, 900)]
    # + strand, 100 upstream/200 downstream
    upstream, downstream, padLeft, padRight = deeptools.heatmapper.chopRegionsFromMiddle(region, left=100, right=200)
    assert(upstream == [(100, 200)])
    assert(downstream == [(300, 400), (800, 900)])
    assert(padLeft == 0)
    assert(padRight == 0)
    # + strand, 250 upstream/300 downstream
    upstream, downstream, padLeft, padRight = deeptools.heatmapper.chopRegionsFromMiddle(region, left=250, right=300)
    assert(upstream == [(0, 200)])
    assert(downstream == [(300, 400), (800, 900)])
    assert(padLeft == 50)
    assert(padRight == 100)
    # - strand, 100 upstream/200 downstream
    upstream, downstream, padLeft, padRight = deeptools.heatmapper.chopRegionsFromMiddle(region, left=200, right=100)
    assert(upstream == [(0, 200)])
    assert(downstream == [(300, 400)])
    assert(padLeft == 0)
    assert(padRight == 0)
    # - strand, 250 upstream/300 downstream
    upstream, downstream, padLeft, padRight = deeptools.heatmapper.chopRegionsFromMiddle(region, left=300, right=250)
    assert(upstream == [(0, 200)])
    assert(downstream == [(300, 400), (800, 900)])
    assert(padLeft == 100)
    assert(padRight == 50)
