import os
import sys
import filecmp
import matplotlib as mpl
import deeptools.plotCoverage
import deeptools.utilities

__author__ = 'Bjoern'

TEST_DATA = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_plotCoverage/"


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

    def test_plotCoverage_default(self):
        args = "--bamfiles {0}test1.bam {0}test2.bam --plotFile /tmp/_test.svg" \
               " --plotFileFormat svg --outRawCounts /tmp/_test.tab".format(TEST_DATA).split()
        deeptools.plotCoverage.main(args)
        assert filecmp.cmp(os.path.join(ROOT, 'outRawCounts_default.tabular'), '/tmp/_test.tab') is True
        if self.run_image_tests and sys.version_info[0] == 2:
            assert self.compare_svg(os.path.join(ROOT, 'plotCoverage_default.svg'), '/tmp/_test.svg')
        os.remove('/tmp/_test.tab')
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
        if res is False:
            os.system("diff {} {}".format(f1, f2))
        os.remove(f1)
        os.remove(f2)
        return res
