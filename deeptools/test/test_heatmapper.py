import os
import sys
import filecmp
import matplotlib as mpl
import deeptools.computeMatrix
import deeptools.plotHeatmap
import deeptools.plotProfile
import deeptools.utilities

__author__ = 'Fidel'

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_heatmapper/"


class TestHeatmapper(object):

    def setUp(self):
        # the tests based on images were done with
        # matplotlib 1.5.0 and will fail if other
        # version is used
        self.run_image_tests = True
        if mpl.__version__ != '1.5.0':
            sys.stderr.write("\nTests based on images are skipped because of "
                             "different matplotlib version ({}) != 1.5.0\n".format(mpl.__version__))
            self.run_image_tests = False

    def test_computeMatrix_reference_point(self):
        args = "reference-point -R {0}/test2.bed -S {0}/test.bw  -b 100 -a 100 " \
               "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1".format(ROOT).split()
        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test.mat.gz')
        assert filecmp.cmp(ROOT + '/master.mat', '/tmp/_test.mat') is True
        os.remove('/tmp/_test.mat')

    def test_computeMatrix_reference_point_missing_data_as_zero(self):
        args = "reference-point -R {0}/test2.bed -S {0}/test.bw  -b 100 -a 100 " \
               "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1 --missingDataAsZero".format(ROOT).split()
        print " ".join(args)
        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test.mat.gz')
        assert filecmp.cmp(ROOT + '/master_nan_to_zero.mat', '/tmp/_test.mat') is True
        os.remove('/tmp/_test.mat')

    def test_computeMatrix_scale_regions(self):
        args = "scale-regions -R {0}/test2.bed -S {0}/test.bw  -b 100 -a 100 -m 100 " \
               "--outFileName /tmp/_test2.mat.gz -bs 10 -p 1".format(ROOT).split()

        deeptools.computeMatrix.main(args)

        os.system('gunzip -f /tmp/_test2.mat.gz')
        assert filecmp.cmp(ROOT + '/master_scale_reg.mat', '/tmp/_test2.mat') is True
        os.remove('/tmp/_test2.mat')

    def test_computeMatrix_multiple_bed(self):
        args = "reference-point -R {0}/group1.bed {0}/group2.bed -S {0}/test.bw  -b 100 -a 100 " \
               "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1".format(ROOT).split()
        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test.mat.gz')
        assert filecmp.cmp(ROOT + '/master_multibed.mat', '/tmp/_test.mat') is True
        os.remove('/tmp/_test.mat')

    def test_computeMatrix_region_extend_over_chr_end(self):
        args = "reference-point -R {0}/group1.bed {0}/group2.bed -S {0}/test.bw  -b 100 -a 500 " \
               "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1".format(ROOT).split()
        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test.mat.gz')
        assert filecmp.cmp(ROOT + '/master_extend_beyond_chr_size.mat', '/tmp/_test.mat') is True
        os.remove('/tmp/_test.mat')

    def test_computeMatrix_unscaled(self):
        args = "scale-regions -S {0}/unscaled.bigWig -R {0}/unscaled.bed -a 300 -b 500 --unscaled5prime 100 --unscaled3prime 50 " \
               "--outFileName /tmp/_test.mat.gz -bs 10 -p 1".format(ROOT).split()
        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test.mat.gz')
        assert filecmp.cmp(ROOT + '/master_unscaled.mat', '/tmp/_test.mat') is True
        os.remove('/tmp/_test.mat')

    def test_plotHeatmap_simple_plot(self):
        """
        Test a simple plot generated using a matrix from
        the following command:

         computeMatrix reference-point -a 100 -b 100 -S {test_path}/test.bw \
          -R {test_path}/test.bed -o /tmp/mat.gz -bs 25

        """
        if self.run_image_tests:
            args = "-m {}/master.mat.gz --outFileName /tmp/_test.svg".format(ROOT).split()
            deeptools.plotHeatmap.main(args)

            # may fail if diff version of matplotlib library is used
            assert self.compare_svg(ROOT + '/master.svg', '/tmp/_test.svg') is True
            os.remove('/tmp/_test.svg')

    def test_plotHeatmap_rename_labels(self):
        if self.run_image_tests:
            args = "-m {}/master.mat.gz --outFileName /tmp/_test2.svg --regionsLabel uno dos".format(ROOT).split()
            deeptools.plotHeatmap.main(args)
            assert self.compare_svg(ROOT + '/master_relabeled.svg', '/tmp/_test2.svg') is True
            os.remove('/tmp/_test2.svg')

    def test_plotHeatmap_scale_regions(self):
        if self.run_image_tests:
            args = "-m {}/master_scale_reg.mat.gz --outFileName /tmp/_test3.svg".format(ROOT).split()
            deeptools.plotHeatmap.main(args)
            assert self.compare_svg(ROOT + '/master_scale_reg.svg', '/tmp/_test3.svg') is True
            os.remove('/tmp/_test3.svg')

    def test_plotHeatmap_multi_bigwig_pergroup(self):
        if self.run_image_tests:
            args = "-m {}/master_multi.mat.gz --perGroup --samplesLabel file1 file2 file3 file4 " \
                   "--outFileName /tmp/_test.svg".format(ROOT).split()
            deeptools.plotHeatmap.main(args)
            assert self.compare_svg(ROOT + '/heatmap_master_multi_pergroup.svg', '/tmp/_test.svg') is True
            os.remove('/tmp/_test.svg')

    def test_plotProfiler(self):
        if self.run_image_tests:
            args = "-m {}/master.mat.gz --outFileName /tmp/_test.svg --regionsLabel uno dos " \
                   "--plotType std".format(ROOT).split()
            deeptools.plotProfile.main(args)
            assert self.compare_svg(ROOT + '/profile_master.svg', '/tmp/_test.svg')
            os.remove('/tmp/_test.svg')

    def test_plotProfiler_heatmap(self):
        if self.run_image_tests:
            args = "-m {}/master.mat.gz --outFileName /tmp/_test.svg --plotType heatmap".format(ROOT).split()
            deeptools.plotProfile.main(args)
            assert self.compare_svg(ROOT + '/profile_master_heatmap.svg', '/tmp/_test.svg')
            os.remove('/tmp/_test.svg')

    def test_plotProfiler_overlapped_lines(self):
        if self.run_image_tests:
            args = "-m {}/master.mat.gz --outFileName /tmp/_test.svg " \
                   "--plotType overlapped_lines --yMin -1".format(ROOT).split()
            deeptools.plotProfile.main(args)
            assert self.compare_svg(ROOT + '/profile_master_overlap_lines.svg', '/tmp/_test.svg')
            os.remove('/tmp/_test.svg')

    def test_plotProfiler_multibigwig(self):
        if self.run_image_tests:
            args = "-m {}/master_multi.mat.gz --outFileName /tmp/_test.svg " \
                   "--numPlotsPerRow 2 --yMax 1.5".format(ROOT).split()
            deeptools.plotProfile.main(args)
            assert self.compare_svg(ROOT + '/profile_master_multi.svg', '/tmp/_test.svg')
            os.remove('/tmp/_test.svg')

    def test_plotProfiler_multibigwig_pergroup(self):
        if self.run_image_tests:
            args = "-m {}/master_multi.mat.gz --outFileName /tmp/_test.svg " \
                   "--perGroup --yMax 1.5".format(ROOT).split()
            deeptools.plotProfile.main(args)
            assert self.compare_svg(ROOT + '/profile_master_multi_pergroup.svg', '/tmp/_test.svg')
            os.remove('/tmp/_test.svg')

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
