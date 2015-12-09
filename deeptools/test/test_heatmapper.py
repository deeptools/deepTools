import os
import subprocess
import filecmp
import deeptools.computeMatrix
import deeptools.plotHeatmap
import deeptools.plotProfile

__author__ = 'Fidel'

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_heatmapper/"


class TestHeatmapper(object):

    def test_computeMatrix_reference_point(self):
        args = "reference-point -R {0}/test2.bed -S {0}/test.bw  -b 100 -a 100 " \
               "--outFileName /tmp/_test.mat.gz  -bs 1 -p 1".format(ROOT).split()
        deeptools.computeMatrix.main(args)
        os.system('gunzip -f /tmp/_test.mat.gz')
        assert filecmp.cmp(ROOT + '/master.mat', '/tmp/_test.mat') is True
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


    def test_plotHeatmap_simple_plot(self):
        """
        Test a simple plot generated using a matrix from
        the following command:

         computeMatrix reference-point -a 100 -b 100 -S {test_path}/test.bw \
          -R {test_path}/test.bed -o /tmp/mat.gz -bs 25

        """
        args = "-m {}/master.mat.gz --outFileName /tmp/_test.svg".format(ROOT).split()
        deeptools.plotHeatmap.main(args)

        # may fail if diff version of matplotlib library is used
        assert self.compare_svg(ROOT + '/master.svg', '/tmp/_test.svg') is True
        os.remove('/tmp/_test.svg')

    def test_plotHeatmap_rename_labels(self):
        args = "-m {}/master.mat.gz --outFileName /tmp/_test2.svg --regionsLabel uno,dos".format(ROOT).split()
        deeptools.plotHeatmap.main(args)
        assert self.compare_svg(ROOT + '/master_relabeled.svg', '/tmp/_test2.svg') is True
        os.remove('/tmp/_test2.svg')

    def test_plotHeatmap_scale_regions(self):
        args = "-m {}/master_scale_reg.mat.gz --outFileName /tmp/_test3.svg".format(ROOT).split()
        deeptools.plotHeatmap.main(args)
        assert self.compare_svg(ROOT + '/master_scale_reg.svg', '/tmp/_test3.svg') is True
        os.remove('/tmp/_test3.svg')

    def test_plotProfiler(self):
        args = "-m {}/master.mat.gz --outFileName /tmp/_test.svg --regionsLabel uno,dos " \
               "--plotType std".format(ROOT).split()
        deeptools.plotProfile.main(args)
        assert self.compare_svg(ROOT + '/profile_master.svg', '/tmp/_test.svg')
        os.remove('/tmp/_test.svg')

        
    @staticmethod
    def compare_svg(file1, file2):
        """
        svg files usually differ on randomly assigned ids.
        This code compares the files ignoring the lines that contain ids

        :return: bool True if files are similar
        """
        try:
            # the diff command is used to compare the files, lines containing the word id are filtered out
            output = subprocess.check_output("/usr/bin/diff  --suppress-common-lines -y "
                                             "{} {} | grep -v id".format(file1, file2), shell=True)
        except subprocess.CalledProcessError as grepexc:
            # if the files are different, diff returns and exit code = 1 that raises
            # the subprocess.CalledProcessError exception
            output = grepexc.output

        if output.strip() == '':
            return True
        else:
            return False
