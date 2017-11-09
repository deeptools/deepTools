import os
import sys
import filecmp
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile

import deeptools.computeMatrix
import deeptools.plotHeatmap
import deeptools.plotProfile
import deeptools.utilities
import json

__author__ = 'Fidel'

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_heatmapper/"

tolerance = 13  # default matplotlib pixed difference tolerance


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


class TestHeatmapper(object):

    def test_computeMatrix_reference_point(self):
        args = "reference-point -R {0}/test2.bed -S {0}/test.bw  -b 100 -a 100 " \
               "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1".format(ROOT).split()
        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test.mat.gz')
        assert cmpMatrices(ROOT + '/master.mat', '/tmp/_test.mat') is True
        os.remove('/tmp/_test.mat')

    def test_computeMatrix_reference_point_center(self):
        args = "reference-point -R {0}/test2.bed -S {0}/test.bw  -b 100 -a 100 --referencePoint center " \
               "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1".format(ROOT).split()
        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test.mat.gz')
        assert cmpMatrices(ROOT + '/master_center.mat', '/tmp/_test.mat') is True
        os.remove('/tmp/_test.mat')

    def test_computeMatrix_reference_point_tes(self):
        args = "reference-point -R {0}/test2.bed -S {0}/test.bw  -b 100 -a 100 --referencePoint TES " \
               "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1".format(ROOT).split()
        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test.mat.gz')
        assert cmpMatrices(ROOT + '/master_TES.mat', '/tmp/_test.mat') is True
        os.remove('/tmp/_test.mat')

    def test_computeMatrix_reference_point_missing_data_as_zero(self):
        args = "reference-point -R {0}/test2.bed -S {0}/test.bw  -b 100 -a 100 " \
               "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1 --missingDataAsZero".format(ROOT).split()
        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test.mat.gz')
        assert cmpMatrices(ROOT + '/master_nan_to_zero.mat', '/tmp/_test.mat') is True
        os.remove('/tmp/_test.mat')

    def test_computeMatrix_scale_regions(self):
        args = "scale-regions -R {0}/test2.bed -S {0}/test.bw  -b 100 -a 100 -m 100 " \
               "--outFileName /tmp/_test2.mat.gz -bs 1 -p 1".format(ROOT).split()

        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test2.mat.gz')
        assert cmpMatrices(ROOT + '/master_scale_reg.mat', '/tmp/_test2.mat') is True
        os.remove('/tmp/_test2.mat')

    def test_computeMatrix_multiple_bed(self):
        args = "reference-point -R {0}/group1.bed {0}/group2.bed -S {0}/test.bw  -b 100 -a 100 " \
               "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1".format(ROOT).split()
        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test.mat.gz')
        assert cmpMatrices(ROOT + '/master_multibed.mat', '/tmp/_test.mat') is True
        os.remove('/tmp/_test.mat')

    def test_computeMatrix_region_extend_over_chr_end(self):
        args = "reference-point -R {0}/group1.bed {0}/group2.bed -S {0}/test.bw  -b 100 -a 500 " \
               "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1".format(ROOT).split()
        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test.mat.gz')
        assert cmpMatrices(ROOT + '/master_extend_beyond_chr_size.mat', '/tmp/_test.mat') is True
        os.remove('/tmp/_test.mat')

    def test_computeMatrix_unscaled(self):
        args = "scale-regions -S {0}/unscaled.bigWig -R {0}/unscaled.bed -a 300 -b 500 --unscaled5prime 100 --unscaled3prime 50 " \
               "--outFileName /tmp/_test.mat.gz -bs 10 -p 1".format(ROOT).split()
        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test.mat.gz')
        assert cmpMatrices(ROOT + '/master_unscaled.mat', '/tmp/_test.mat') is True
        os.remove('/tmp/_test.mat')

    def test_computeMatrix_gtf(self):
        args = "scale-regions -S {0}../test_data/test1.bw.bw -R {0}../test_data/test.gtf -a 300 -b 500 --unscaled5prime 20 --unscaled3prime 50 " \
               "--outFileName /tmp/_test_gtf.mat.gz -bs 10 -p 1".format(ROOT).split()
        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test_gtf.mat.gz')
        assert cmpMatrices(ROOT + '/master_gtf.mat', '/tmp/_test_gtf.mat') is True
        os.remove('/tmp/_test_gtf.mat')

    def test_computeMatrix_metagene(self):
        args = "scale-regions -S {0}../test_data/test1.bw.bw -R {0}../test_data/test.gtf -a 300 -b 500 --unscaled5prime 20 --unscaled3prime 50 " \
               "--outFileName /tmp/_test_metagene.mat.gz -bs 10 -p 1 --metagene".format(ROOT).split()
        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test_metagene.mat.gz')
        assert cmpMatrices(ROOT + '/master_metagene.mat', '/tmp/_test_metagene.mat') is True
        os.remove('/tmp/_test_metagene.mat')

    def test_plotHeatmap_simple_plot(self):
        """
        Test a simple plot generated using a matrix from
        the following command:

         computeMatrix reference-point -a 100 -b 100 -S {test_path}/test.bw \
          -R {test_path}/test.bed -o /tmp/mat.gz -bs 25

        """
        outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
        args = "-m {}/master.mat.gz --outFileName {}".format(ROOT, outfile.name).split()
        deeptools.plotHeatmap.main(args)
        res = compare_images(ROOT + '/master.png', outfile.name, tolerance)
        assert res is None, res
        os.remove(outfile.name)

    def test_plotHeatmap_rename_labels(self):
        outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)

        args = "-m {}/master.mat.gz --outFileName {} --regionsLabel uno dos".format(ROOT, outfile.name).split()
        deeptools.plotHeatmap.main(args)
        res = compare_images(ROOT + '/master_relabeled.png', outfile.name, tolerance)
        assert res is None, res
        os.remove(outfile.name)

    def test_plotHeatmap_scale_regions(self):
        outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
        args = "-m {}/master_scale_reg.mat.gz --outFileName {}".format(ROOT, outfile.name).split()
        deeptools.plotHeatmap.main(args)
        res = compare_images(ROOT + '/master_scale_reg.png', outfile.name, tolerance)
        assert res is None, res
        os.remove(outfile.name)

    def test_plotHeatmap_multi_bigwig_pergroup(self):
        outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
        args = "-m {}/master_multi.mat.gz --perGroup --samplesLabel file1 file2 file3 file4 " \
               "--outFileName {}".format(ROOT, outfile.name).split()
        deeptools.plotHeatmap.main(args)
        res = compare_images(ROOT + '/heatmap_master_multi_pergroup.png', outfile.name, tolerance)
        assert res is None, res
        os.remove(outfile.name)

    def test_plotHeatmap_multiple_colors_muti_scales(self):
        outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
        args = "-m {}/master_multi.mat.gz --colorList white,blue white,red --zMin 1 0 --zMax 4 5 " \
               "--outFileName {}".format(ROOT, outfile.name).split()
        deeptools.plotHeatmap.main(args)
        res = compare_images(ROOT + '/heatmap_master_multi_color.png', outfile.name, tolerance)
        assert res is None, res
        os.remove(outfile.name)

    def test_plotHeatmap_multiple_colormap_no_boxes(self):
        outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
        args = "-m {}/master_multi.mat.gz --colorMap Reds binary terrain --boxAroundHeatmaps no " \
               "--outFileName {}".format(ROOT, outfile.name).split()
        deeptools.plotHeatmap.main(args)
        res = compare_images(ROOT + '/heatmap_master_multi_colormap_no_box.png', outfile.name, tolerance)
        assert res is None, res
        os.remove(outfile.name)

    def test_plotHeatmap_interpolation(self):
        outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
        args = "-m {}/large_matrix.mat.gz --interpolation bilinear " \
               "--outFileName {}".format(ROOT, outfile.name).split()
        deeptools.plotHeatmap.main(args)
        print " ".join(args)
        res = compare_images(ROOT + '/heatmap_master_interpolation_bilinear.png', outfile.name, tolerance)
        assert res is None, res
        os.remove(outfile.name)

    def test_plotProfiler(self):
        outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
        args = "-m {}/master.mat.gz --outFileName {} --regionsLabel uno dos " \
               "--plotType std".format(ROOT, outfile.name).split()
        deeptools.plotProfile.main(args)
        res = compare_images(ROOT + '/profile_master.png', outfile.name, tolerance)
        assert res is None, res
        os.remove(outfile.name)

    def test_plotProfiler_heatmap(self):
        outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
        args = "-m {}/master.mat.gz --outFileName {} --plotType heatmap".format(ROOT, outfile.name).split()
        deeptools.plotProfile.main(args)
        res = compare_images(ROOT + '/profile_master_heatmap.png', outfile.name, tolerance)
        assert res is None, res
        os.remove(outfile.name)

    def test_plotProfiler_overlapped_lines(self):
        outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
        args = "-m {}/master.mat.gz --outFileName {} " \
               "--plotType overlapped_lines --yMin -1".format(ROOT, outfile.name).split()
        deeptools.plotProfile.main(args)
        res = compare_images(ROOT + '/profile_master_overlap_lines.png', outfile.name, tolerance)
        assert res is None, res
        os.remove(outfile.name)

    def test_plotProfiler_multibigwig(self):
        outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
        args = "-m {}/master_multi.mat.gz --outFileName {} " \
               "--numPlotsPerRow 2 --yMax 1.5".format(ROOT, outfile.name).split()
        deeptools.plotProfile.main(args)
        res = compare_images(ROOT + '/profile_master_multi.png', outfile.name, tolerance)
        assert res is None, res
        os.remove(outfile.name)

    def test_plotProfiler_multibigwig_pergroup(self):
        outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
        args = "-m {}/master_multi.mat.gz --outFileName {} " \
               "--perGroup --yMax 1.5".format(ROOT, outfile.name).split()
        deeptools.plotProfile.main(args)
        res = compare_images(ROOT + '/profile_master_multi_pergroup.png', outfile.name, tolerance)
        assert res is None, res
        os.remove(outfile.name)

    def test_chopRegions_body(self):
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

    def test_chopRegions_TSS(self):
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

    def test_chopRegions_TES(self):
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

    def test_chopRegionsFromMiddle(self):
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

    @staticmethod
    def compare_svg(file1, file2):
        """
        svg files usually differ on randomly assigned ids and xlink:href tags
        This code compares the files ignoring the lines that contain ids

        :return: bool True if files are similar
        """
        f1 = deeptools.utilities.getTempFileName(suffix='.svg')
        f2 = deeptools.utilities.getTempFileName(suffix='.svg')
        # remove xlink:href, id and url attributes
        os.system('cat {} | perl -lane \'s/xlink:href=".+?"//g; s/id=".+?"//g; s/"url\(.+?\)"//g; print $_\' > {}'.format(file1, f1))
        os.system('cat {} | perl -lane \'s/xlink:href=".+?"//g; s/id=".+?"//g; s/"url\(.+?\)"//g; print $_\' > {}'.format(file2, f2))
        res = filecmp.cmp(f1, f2)
        os.remove(f1)
        os.remove(f2)
        return res
