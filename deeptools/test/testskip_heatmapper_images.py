import os
import matplotlib
matplotlib.use('Agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile

import deeptools.computeMatrix
import deeptools.plotHeatmap
import deeptools.plotProfile
import deeptools.utilities

__author__ = 'Fidel'

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_heatmapper/"
tolerance = 30

skip = False
if matplotlib.__version__ != "3.1.1":
    skip = True


def test_plotHeatmap_simple_plot():
    """
    Test a simple plot generated using a matrix from
    the following command:

     computeMatrix reference-point -a 100 -b 100 -S {test_path}/test.bw \
      -R {test_path}/test.bed -o /tmp/mat.gz -bs 25

    """
    if skip:
        return
    outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
    args = "-m {}/master.mat.gz --outFileName {}".format(ROOT, outfile.name).split()
    deeptools.plotHeatmap.main(args)
    res = compare_images(ROOT + '/master.png', outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_plotHeatmap_rename_labels():
    if skip:
        return
    outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)

    args = "-m {}/master.mat.gz --outFileName {} --regionsLabel uno dos".format(ROOT, outfile.name).split()
    deeptools.plotHeatmap.main(args)
    res = compare_images(ROOT + '/master_relabeled.png', outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_plotHeatmap_scale_regions():
    if skip:
        return
    outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
    args = "-m {}/master_scale_reg.mat.gz --outFileName {}".format(ROOT, outfile.name).split()
    deeptools.plotHeatmap.main(args)
    res = compare_images(ROOT + '/master_scale_reg.png', outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_plotHeatmap_multi_bigwig_pergroup():
    if skip:
        return
    outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
    args = "-m {}/master_multi.mat.gz --perGroup --samplesLabel file1 file2 file3 file4 " \
           "--outFileName {}".format(ROOT, outfile.name).split()
    deeptools.plotHeatmap.main(args)
    res = compare_images(ROOT + '/heatmap_master_multi_pergroup.png', outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_plotHeatmap_multiple_colors_muti_scales():
    if skip:
        return
    outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
    args = "-m {}/master_multi.mat.gz --colorList white,blue white,red --zMin 1 0 --zMax 4 5 " \
           "--outFileName {}".format(ROOT, outfile.name).split()
    deeptools.plotHeatmap.main(args)
    res = compare_images(ROOT + '/heatmap_master_multi_color.png', outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_plotHeatmap_multiple_colormap_no_boxes():
    if skip:
        return
    outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
    args = "-m {}/master_multi.mat.gz --colorMap Reds binary terrain --boxAroundHeatmaps no " \
           "--outFileName {}".format(ROOT, outfile.name).split()
    deeptools.plotHeatmap.main(args)
    res = compare_images(ROOT + '/heatmap_master_multi_colormap_no_box.png', outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_plotHeatmap_interpolation():
    if skip:
        return
    outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
    args = "-m {}/large_matrix.mat.gz --interpolation bilinear " \
           "--outFileName {}".format(ROOT, outfile.name).split()
    deeptools.plotHeatmap.main(args)
    res = compare_images(ROOT + '/heatmap_master_interpolation_bilinear.png', outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_plotProfiler():
    if skip:
        return
    outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
    args = "-m {}/master.mat.gz --outFileName {} --regionsLabel uno dos " \
           "--plotType std".format(ROOT, outfile.name).split()
    deeptools.plotProfile.main(args)
    res = compare_images(ROOT + '/profile_master.png', outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_plotProfiler_heatmap():
    if skip:
        return
    outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
    args = "-m {}/master.mat.gz --outFileName {} --plotType heatmap".format(ROOT, outfile.name).split()
    deeptools.plotProfile.main(args)
    res = compare_images(ROOT + '/profile_master_heatmap.png', outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_plotProfiler_overlapped_lines():
    if skip:
        return
    outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
    args = "-m {}/master.mat.gz --outFileName {} " \
           "--plotType overlapped_lines --yMin -1".format(ROOT, outfile.name).split()
    deeptools.plotProfile.main(args)
    res = compare_images(ROOT + '/profile_master_overlap_lines.png', outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_plotProfiler_multibigwig():
    if skip:
        return
    outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
    args = "-m {}/master_multi.mat.gz --outFileName {} " \
           "--numPlotsPerRow 2 --yMax 1.5".format(ROOT, outfile.name).split()
    deeptools.plotProfile.main(args)
    res = compare_images(ROOT + '/profile_master_multi.png', outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_plotProfiler_multibigwig_pergroup():
    if skip:
        return
    outfile = NamedTemporaryFile(suffix='.png', prefix='plotHeatmap_test_', delete=False)
    args = "-m {}/master_multi.mat.gz --outFileName {} " \
           "--perGroup --yMax 1.5".format(ROOT, outfile.name).split()
    deeptools.plotProfile.main(args)
    res = compare_images(ROOT + '/profile_master_multi_pergroup.png', outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)
