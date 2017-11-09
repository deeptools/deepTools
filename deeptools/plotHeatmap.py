#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division

import argparse
from collections import OrderedDict
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.gridspec as gridspec
from matplotlib import ticker

import sys

# own modules
from deeptools import parserCommon
from deeptools import heatmapper
from deeptools.heatmapper_utilities import plot_single, getProfileTicks

debug = 0
old_settings = np.seterr(all='ignore')
plt.ioff()


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        parents=[parserCommon.heatmapperMatrixArgs(),
                 parserCommon.heatmapperOutputArgs(mode='heatmap'),
                 parserCommon.heatmapperOptionalArgs(mode='heatmap')],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='This tool creates a heatmap for '
        'scores associated with genomic regions. '
        'The program requires a matrix file '
        'generated by the tool ``computeMatrix``.',
        epilog='An example usage is: plotHeatmap -m <matrix file>',
        add_help=False)

    return parser


def process_args(args=None):
    args = parse_arguments().parse_args(args)

    args.heatmapHeight = args.heatmapHeight if args.heatmapHeight > 3 and args.heatmapHeight <= 100 else 10

    if not matplotlib.colors.is_color_like(args.missingDataColor):
        exit("The value {0}  for --missingDataColor is not valid".format(args.missingDataColor))

    args.boxAroundHeatmaps = True if args.boxAroundHeatmaps == 'yes' else False

    return args


def get_heatmap_ticks(hm, reference_point_label, startLabel, endLabel):
    """
    returns the position and labelling of the xticks that
    correspond to the heatmap
    """
    w = hm.parameters['bin size']
    b = hm.parameters['upstream']
    a = hm.parameters['downstream']
    c = hm.parameters.get('unscaled 5 prime', 0)
    d = hm.parameters.get('unscaled 3 prime', 0)
    m = hm.parameters['body']

    if b < 1e5:
        quotient = 1000
        symbol = 'Kb'
    else:
        quotient = 1e6
        symbol = 'Mb'

    if m == 0:
        xticks = [(k / w) for k in [0, b, b + a]]
        xticks_label = ['{0:.1f}'.format(-(float(b) / quotient)), reference_point_label,
                        '{0:.1f}{1}'.format(float(a) / quotient, symbol)]

    else:
        xticks_values = [0]
        xticks_label = []

        # only if upstream region is set, add a x tick
        if b > 0:
            xticks_values.append(b)
            xticks_label.append('{0:.1f}'.format(-(float(b) / quotient)))

        xticks_label.append(startLabel)

        # 5 prime unscaled tick/label, if needed
        if c > 0:
            xticks_values.append(b + c)
            xticks_label.append('')

        # set the x tick for the body parameter, regardless if
        # upstream is 0 (not set)
        if d > 0:
            xticks_values.append(b + c + m)
            xticks_label.append('')

        xticks_values.append(b + c + m + d)
        xticks_label.append(endLabel)

        if a > 0:
            xticks_values.append(b + c + m + d + a)
            xticks_label.append('{0:.1f}{1}'.format(float(a) / quotient, symbol))

        xticks = [k / w for k in xticks_values]

    return xticks, xticks_label


def prepare_layout(hm_matrix, heatmapsize, showSummaryPlot, showColorbar, perGroup, colorbar_position):
    """
    prepare the plot layout
    as a grid having as many rows
    as samples (+1 for colobar)
    and as many rows as groups (or clusters) (+1 for profile plot)
    """
    heatmapwidth, heatmapheight = heatmapsize

    numcols = hm_matrix.get_num_samples()
    numrows = hm_matrix.get_num_groups()
    if perGroup:
        numcols, numrows = numrows, numcols

    # the rows have different size depending
    # on the number of regions contained in the
    if perGroup:
        # heatmap
        height_ratio = np.array([np.amax(np.diff(hm_matrix.group_boundaries))] * numrows)
        # scale ratio to sum = heatmapheight
        height_ratio = heatmapheight * (height_ratio.astype(float) / height_ratio.sum())
    else:
        # heatmap
        height_ratio = np.diff(hm_matrix.group_boundaries)
        # scale ratio to sum = heatmapheight
        height_ratio = heatmapheight * (height_ratio.astype(float) / height_ratio.sum())

    # convert the height_ratio from numpy array back to list
    height_ratio = height_ratio.tolist()
    # the width ratio is equal for all heatmaps
    width_ratio = [heatmapwidth] * numcols

    if showColorbar:
        if colorbar_position == 'below':
            numrows += 2  # a spacer needs to be added to avoid overlaps
            height_ratio += [4 / 2.54]  # spacer
            height_ratio += [1 / 2.54]
        else:
            numcols += 1
            width_ratio += [1 / 2.54]

    if showSummaryPlot:
        numrows += 2  # plus 2 because a spacer is added
        # make height of summary plot
        # proportional to the width of heatmap
        sumplot_height = heatmapwidth
        spacer_height = heatmapwidth / 8
        # scale height_ratios to convert from row
        # numbers to heatmapheigt fractions
        height_ratio = np.concatenate([[sumplot_height, spacer_height], height_ratio])

    grids = gridspec.GridSpec(numrows, numcols, height_ratios=height_ratio, width_ratios=width_ratio)

    return grids


def addProfilePlot(hm, plt, fig, grids, iterNum, iterNum2, perGroup, averageType, xticks, xtickslabel, yAxisLabel, color_list, yMin, yMax, wspace, hspace, colorbar_position):
    """
    A function to add profile plots to the given figure, possibly in a custom grid subplot which mimics a tight layout (if wspace and hspace are not None)
    """
    if wspace is not None and hspace is not None:
        if colorbar_position == 'side':
            gridsSub = gridspec.GridSpecFromSubplotSpec(1, iterNum, subplot_spec=grids[0, :-1], wspace=wspace, hspace=hspace)
        else:
            gridsSub = gridspec.GridSpecFromSubplotSpec(1, iterNum, subplot_spec=grids[0, :], wspace=wspace, hspace=hspace)

    ax_list = []
    for sample_id in range(iterNum):
        if perGroup:
            title = hm.matrix.group_labels[sample_id]
        else:
            title = hm.matrix.sample_labels[sample_id]
        if sample_id > 0 and len(yMin) == 1 and len(yMax) == 1:
            ax_profile = fig.add_subplot(grids[0, sample_id], sharey=ax_list[0])
        else:
            if wspace is not None and hspace is not None:
                ax_profile = fig.add_subplot(gridsSub[0, sample_id])
            else:
                ax_profile = fig.add_subplot(grids[0, sample_id])

        ax_profile.set_title(title)
        for group in range(iterNum2):
            if perGroup:
                sub_matrix = hm.matrix.get_matrix(sample_id, group)
                line_label = sub_matrix['sample']
            else:
                sub_matrix = hm.matrix.get_matrix(group, sample_id)
                line_label = sub_matrix['group']
            plot_single(ax_profile, sub_matrix['matrix'],
                        averageType,
                        color_list[group],
                        line_label,
                        plot_type='simple')

        if sample_id > 0 and len(yMin) == 1 and len(yMax) == 1:
            plt.setp(ax_profile.get_yticklabels(), visible=False)

        if sample_id == 0 and yAxisLabel != '':
            ax_profile.set_ylabel(yAxisLabel)
        if np.ceil(max(xticks)) != float(sub_matrix['matrix'].shape[1]):
            tickscale = float(sub_matrix['matrix'].shape[1]) / max(xticks)
            xticks_use = [x * tickscale for x in xticks]
            ax_profile.axes.set_xticks(xticks_use)
        else:
            ax_profile.axes.set_xticks(xticks)
        ax_profile.axes.set_xticklabels(xtickslabel)
        ax_list.append(ax_profile)

        # align the first and last label
        # such that they don't fall off
        # the heatmap sides
        ticks = ax_profile.xaxis.get_major_ticks()
        ticks[0].label1.set_horizontalalignment('left')
        ticks[-1].label1.set_horizontalalignment('right')

        # It turns out that set_ylim only takes np.float64s
        localYMin = yMin[sample_id % len(yMin)]
        localYMax = yMax[sample_id % len(yMax)]
        lims = ax_list[0].get_ylim()
        if localYMin:
            lims = (np.float64(localYMin), lims[1])
        if localYMax:
            lims = (lims[0], np.float64(localYMax))
        if lims[0] >= lims[1]:
            lims = (lims[0], lims[0] + 1)
        if len(yMin) == 1 and len(yMax) == 1:
            ax_list[0].set_ylim(lims)
        else:
            ax_list[-1].set_ylim(lims)
    return ax_list


def plotMatrix(hm, outFileName,
               colorMapDict={'colorMap': ['binary'], 'missingDataColor': 'black', 'alpha': 1.0},
               plotTitle='',
               xAxisLabel='', yAxisLabel='', regionsLabel='',
               zMin=None, zMax=None,
               yMin=None, yMax=None,
               averageType='median',
               reference_point_label='TSS',
               startLabel='TSS', endLabel="TES",
               heatmapHeight=25,
               heatmapWidth=7.5,
               perGroup=False, whatToShow='plot, heatmap and colorbar',
               image_format=None,
               legend_location='upper-left',
               box_around_heatmaps=True,
               dpi=200,
               interpolation_method='auto'):

    matrix_flatten = None
    if zMin is None:
        matrix_flatten = hm.matrix.flatten()
        # try to avoid outliers by using np.percentile
        zMin = np.percentile(matrix_flatten, 1.0)
        if np.isnan(zMin):
            zMin = [None]
        else:
            zMin = [zMin]  # convert to list to support multiple entries

    if zMax is None:
        if matrix_flatten is None:
            matrix_flatten = hm.matrix.flatten()
        # try to avoid outliers by using np.percentile
        zMax = np.percentile(matrix_flatten, 98.0)
        if np.isnan(zMax) or zMax <= zMin[0]:
            zMax = [None]
        else:
            zMax = [zMax]

    if yMin is None:
        yMin = [None]
    if yMax is None:
        yMax = [None]
    if not isinstance(yMin, list):
        yMin = [yMin]
    if not isinstance(yMax, list):
        yMax = [yMax]

    plt.rcParams['font.size'] = 8.0
    fontP = FontProperties()

    showSummaryPlot = False
    showColorbar = False

    if whatToShow == 'plot and heatmap':
        showSummaryPlot = True
    elif whatToShow == 'heatmap and colorbar':
        showColorbar = True
    elif whatToShow == 'plot, heatmap and colorbar':
        showSummaryPlot = True
        showColorbar = True

    # colormap for the heatmap
    if colorMapDict['colorMap']:
        cmap = []
        for color_map in colorMapDict['colorMap']:
            cmap.append(plt.get_cmap(color_map))
            cmap[-1].set_bad(colorMapDict['missingDataColor'])  # nans are printed using this color

    if colorMapDict['colorList'] and len(colorMapDict['colorList']) > 0:
        # make a cmap for each color list given
        cmap = []
        for color_list in colorMapDict['colorList']:
            cmap.append(matplotlib.colors.LinearSegmentedColormap.from_list(
                'my_cmap', color_list.replace(' ', '').split(","), N=colorMapDict['colorNumber']))
            cmap[-1].set_bad(colorMapDict['missingDataColor'])  # nans are printed using this color

    if len(cmap) > 1 or len(zMin) > 1 or len(zMax) > 1:
        # position color bar below heatmap when more than one
        # heatmap color is given
        colorbar_position = 'below'
    else:
        colorbar_position = 'side'

    grids = prepare_layout(hm.matrix, (heatmapWidth, heatmapHeight),
                           showSummaryPlot, showColorbar, perGroup, colorbar_position)

    # figsize: w,h tuple in inches
    figwidth = heatmapWidth / 2.54
    figheight = heatmapHeight / 2.54
    if showSummaryPlot:
        # the summary plot ocupies a height
        # equal to the fig width
        figheight += figwidth

    numsamples = hm.matrix.get_num_samples()
    if perGroup:
        num_cols = hm.matrix.get_num_groups()
    else:
        num_cols = numsamples
    total_figwidth = figwidth * num_cols
    if showColorbar:
        if colorbar_position == 'below':
            figheight += 1 / 2.54
        else:
            total_figwidth += 1 / 2.54

    fig = plt.figure(figsize=(total_figwidth, figheight))

    xticks, xtickslabel = getProfileTicks(hm, reference_point_label, startLabel, endLabel)

    xticks_heat, xtickslabel_heat = get_heatmap_ticks(hm, reference_point_label, startLabel, endLabel)
    fig.suptitle(plotTitle, y=1 - (0.06 / figheight))

    # color map for the summary plot (profile) on top of the heatmap
    cmap_plot = plt.get_cmap('jet')
    numgroups = hm.matrix.get_num_groups()
    if perGroup:
        color_list = cmap_plot(np.arange(hm.matrix.get_num_samples()) / hm.matrix.get_num_samples())
    else:
        color_list = cmap_plot(np.arange(numgroups) / numgroups)
    alpha = colorMapDict['alpha']

    # check if matrix is reference-point based using the upstream >0 value
    # and is sorted by region length. If this is
    # the case, prepare the data to plot a border at the regions end
    if hm.parameters['upstream'] > 0 and \
            hm.matrix.sort_using == 'region_length' and \
            hm.matrix.sort_method != 'no':

            _regions = hm.matrix.get_regions()
            regions_length_in_bins = []
            for _group in _regions:
                _reg_len = []
                for ind_reg in _group:
                    if isinstance(ind_reg, dict):
                        _len = ind_reg['end'] - ind_reg['start']
                    else:
                        _len = sum([x[1] - x[0] for x in ind_reg[1]])
                    _reg_len.append((hm.parameters['upstream'] + _len) / hm.parameters['bin size'])
                regions_length_in_bins.append(_reg_len)
    else:
        regions_length_in_bins = None

    # plot the profiles on top of the heatmaps
    if showSummaryPlot:
        if perGroup:
            iterNum = numgroups
            iterNum2 = hm.matrix.get_num_samples()
        else:
            iterNum = hm.matrix.get_num_samples()
            iterNum2 = numgroups
        ax_list = addProfilePlot(hm, plt, fig, grids, iterNum, iterNum2, perGroup, averageType, xticks, xtickslabel, yAxisLabel, color_list, yMin, yMax, None, None, colorbar_position)
        if len(yMin) > 1 or len(yMax) > 1:
            # replot with a tight layout
            import matplotlib.tight_layout as tl
            specList = tl.get_subplotspec_list(fig.axes, grid_spec=grids)
            renderer = tl.get_renderer(fig)
            kwargs = tl.get_tight_layout_figure(fig, fig.axes, specList, renderer, pad=1.08)

            for ax in ax_list:
                fig.delaxes(ax)

            ax_list = addProfilePlot(hm, plt, fig, grids, iterNum, iterNum2, perGroup, averageType, xticks, xtickslabel, yAxisLabel, color_list, yMin, yMax, kwargs['wspace'], kwargs['hspace'], colorbar_position)

        if legend_location != 'none':
            ax_list[-1].legend(loc=legend_location.replace('-', ' '), ncol=1, prop=fontP,
                               frameon=False, markerscale=0.5)

    first_group = 0  # helper variable to place the title per sample/group
    for sample in range(hm.matrix.get_num_samples()):
        sample_idx = sample
        for group in range(numgroups):
            group_idx = group
            # add the respective profile to the
            # summary plot
            sub_matrix = hm.matrix.get_matrix(group, sample)
            if showSummaryPlot:
                if perGroup:
                    sample_idx = sample + 2  # plot + spacer
                else:
                    group += 2  # plot + spacer
                first_group = 1

            if perGroup:
                ax = fig.add_subplot(grids[sample_idx, group])
                # the remainder (%) is used to iterate
                # over the available color maps (cmap).
                # if the user only provided, lets say two
                # and there are 10 groups, colormaps they are reused every
                # two groups.
                cmap_idx = group_idx % len(cmap)
                zmin_idx = group_idx % len(zMin)
                zmax_idx = group_idx % len(zMax)
            else:
                ax = fig.add_subplot(grids[group, sample])
                # see above for the use of '%'
                cmap_idx = sample % len(cmap)
                zmin_idx = sample % len(zMin)
                zmax_idx = sample % len(zMax)

            if group == first_group and not showSummaryPlot and not perGroup:
                title = hm.matrix.sample_labels[sample]
                ax.set_title(title)

            if box_around_heatmaps is False:
                # Turn off the boxes around the individual heatmaps
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.spines['left'].set_visible(False)
            rows, cols = sub_matrix['matrix'].shape
            # if the number of rows is too large, then the 'nearest' method simply
            # drops rows. A better solution is to relate the threshold to the DPI of the image
            if interpolation_method == 'auto':
                if rows >= 1000:
                    interpolation_method = 'bilinear'
                else:
                    interpolation_method = 'nearest'

            # if np.clip is not used, then values of the matrix that exceed the zmax limit are
            # highlighted. Usually, a significant amount of pixels are equal or above the zmax and
            # the default behaviour produces images full of large highlighted dots.
            # If interpolation='nearest' is used, this has no effect
            sub_matrix['matrix'] = np.clip(sub_matrix['matrix'], zMin[zmin_idx], zMax[zmax_idx])
            img = ax.imshow(sub_matrix['matrix'],
                            aspect='auto',
                            interpolation=interpolation_method,
                            origin='upper',
                            vmin=zMin[zmin_idx],
                            vmax=zMax[zmax_idx],
                            cmap=cmap[cmap_idx],
                            alpha=alpha,
                            extent=[0, cols, rows, 0])
            img.set_rasterized(True)
            # plot border at the end of the regions
            # if ordered by length
            if regions_length_in_bins is not None:
                x_lim = ax.get_xlim()
                y_lim = ax.get_ylim()

                ax.plot(regions_length_in_bins[group_idx],
                        np.arange(len(regions_length_in_bins[group_idx])),
                        '--', color='black', linewidth=0.5, dashes=(3, 2))
                ax.set_xlim(x_lim)
                ax.set_ylim(y_lim)

            if perGroup:
                ax.axes.set_xlabel(sub_matrix['group'])
                if sample < hm.matrix.get_num_samples() - 1:
                    ax.axes.get_xaxis().set_visible(False)
            else:
                ax.axes.get_xaxis().set_visible(False)
                ax.axes.set_xlabel(xAxisLabel)
            ax.axes.set_yticks([])
            if perGroup and group == 0:
                ax.axes.set_ylabel(sub_matrix['sample'])
            elif not perGroup and sample == 0:
                ax.axes.set_ylabel(sub_matrix['group'])

            # add labels to last block in a column
            if (perGroup and sample == numsamples - 1) or \
               (not perGroup and group_idx == numgroups - 1):

                # add xticks to the bottom heatmap (last group)
                ax.axes.get_xaxis().set_visible(True)
                if np.ceil(max(xticks_heat)) != float(sub_matrix['matrix'].shape[1]):
                    tickscale = float(sub_matrix['matrix'].shape[1]) / max(xticks_heat)
                    xticks_heat_use = [x * tickscale for x in xticks_heat]
                    ax.axes.set_xticks(xticks_heat_use)
                else:
                    ax.axes.set_xticks(xticks_heat)
                ax.axes.set_xticklabels(xtickslabel_heat, size=8)

                # align the first and last label
                # such that they don't fall off
                # the heatmap sides
                ticks = ax.xaxis.get_major_ticks()
                ticks[0].label1.set_horizontalalignment('left')
                ticks[-1].label1.set_horizontalalignment('right')

                ax.get_xaxis().set_tick_params(
                    which='both',
                    top='off',
                    direction='out')

                if showColorbar and colorbar_position == 'below':
                    # draw a colormap per each heatmap below the last block
                    if perGroup:
                        col = group_idx
                    else:
                        col = sample
                    ax = fig.add_subplot(grids[-1, col])
                    tick_locator = ticker.MaxNLocator(nbins=3)
                    cbar = fig.colorbar(img, cax=ax, alpha=alpha, orientation='horizontal', ticks=tick_locator)
                    labels = cbar.ax.get_xticklabels()
                    ticks = cbar.ax.get_xticks()
                    if ticks[0] == 0:
                        # if the label is at the start of the colobar
                        # move it a bit inside to avoid overlapping
                        # with other labels
                        labels[0].set_horizontalalignment('left')
                    if ticks[-1] == 1:
                        # if the label is at the end of the colobar
                        # move it a bit inside to avoid overlapping
                        # with other labels
                        labels[-1].set_horizontalalignment('right')
                    # cbar.ax.set_xticklabels(labels, rotation=90)

    if showColorbar and colorbar_position != 'below':
        if showSummaryPlot:
            # we don't want to colorbar to extend
            # over the profiles and spacer top rows
            grid_start = 2
        else:
            grid_start = 0

        ax = fig.add_subplot(grids[grid_start:, -1])
        fig.colorbar(img, cax=ax, alpha=alpha)

    if box_around_heatmaps:
        plt.subplots_adjust(wspace=0.10, hspace=0.025, top=0.85, bottom=0, left=0.04, right=0.96)
    else:
        #  When no box is plotted the space between heatmaps is reduced
        plt.subplots_adjust(wspace=0.05, hspace=0.01, top=0.85, bottom=0, left=0.04, right=0.96)

    plt.savefig(outFileName, bbox_inches='tight', pdd_inches=0, dpi=dpi, format=image_format)
    plt.close()


def mergeSmallGroups(matrixDict):
    group_lengths = [len(x) for x in matrixDict.values()]
    min_group_length = sum(group_lengths) * 0.01

    to_merge = []
    i = 0
    _mergedHeatMapDict = OrderedDict()

    for label, ma in matrixDict.items():
        # merge small groups together
        # otherwise visualization is impaired
        if group_lengths[i] > min_group_length:
            if len(to_merge):
                to_merge.append(label)
                new_label = " ".join(to_merge)
                new_ma = np.concatenate([matrixDict[item]
                                        for item in to_merge], axis=0)
            else:
                new_label = label
                new_ma = matrixDict[label]

            _mergedHeatMapDict[new_label] = new_ma
            to_merge = []
        else:
            to_merge.append(label)
        i += 1
    if len(to_merge) > 1:
        new_label = " ".join(to_merge)
        new_ma = np.array()
        for item in to_merge:
            new_ma = np.concatenate([new_ma, matrixDict[item]])
        _mergedHeatMapDict[new_label] = new_ma

    return _mergedHeatMapDict


def main(args=None):
    args = process_args(args)
    hm = heatmapper.heatmapper()
    matrix_file = args.matrixFile.name
    args.matrixFile.close()
    hm.read_matrix_file(matrix_file)

    if args.kmeans is not None:
        hm.matrix.hmcluster(args.kmeans, method='kmeans')
    else:
        if args.hclust is not None:
            print("Performing hierarchical clustering."
                  "Please note that it might be very slow for large datasets.\n")
            hm.matrix.hmcluster(args.hclust, method='hierarchical')

    group_len_ratio = np.diff(hm.matrix.group_boundaries) / len(hm.matrix.regions)
    if np.any(group_len_ratio < 5.0 / 1000):
        problem = np.flatnonzero(group_len_ratio < 5.0 / 1000)
        sys.stderr.write("WARNING: Group '{}' is too small for plotting, you might want to remove it. "
                         "There will likely be an error message from matplotlib regarding this "
                         "below.\n".format(hm.matrix.group_labels[problem[0]]))

    if args.regionsLabel:
        hm.matrix.set_group_labels(args.regionsLabel)

    if args.samplesLabel and len(args.samplesLabel):
        hm.matrix.set_sample_labels(args.samplesLabel)

    if args.sortRegions != 'no':
        sortUsingSamples = []
        if args.sortUsingSamples is not None:
            for i in args.sortUsingSamples:
                if (i > 0 and i <= hm.matrix.get_num_samples()):
                    sortUsingSamples.append(i - 1)
                else:
                    exit("The value {0} for --sortSamples is not valid. Only values from 1 to {1} are allowed.".format(args.sortUsingSamples, hm.matrix.get_num_samples()))
            print('Samples used for ordering within each group: ', sortUsingSamples)

        hm.matrix.sort_groups(sort_using=args.sortUsing,
                              sort_method=args.sortRegions,
                              sample_list=sortUsingSamples)

    if args.outFileNameMatrix:
        hm.save_matrix_values(args.outFileNameMatrix)

    if args.outFileSortedRegions:
        hm.save_BED(args.outFileSortedRegions)

    colormap_dict = {'colorMap': args.colorMap,
                     'colorList': args.colorList,
                     'colorNumber': args.colorNumber,
                     'missingDataColor': args.missingDataColor,
                     'alpha': args.alpha}

    plotMatrix(hm,
               args.outFileName,
               colormap_dict, args.plotTitle,
               args.xAxisLabel, args.yAxisLabel, args.regionsLabel,
               args.zMin, args.zMax,
               args.yMin, args.yMax,
               args.averageTypeSummaryPlot,
               args.refPointLabel,
               args.startLabel,
               args.endLabel,
               args.heatmapHeight,
               args.heatmapWidth,
               args.perGroup,
               args.whatToShow,
               image_format=args.plotFileFormat,
               legend_location=args.legendLocation,
               box_around_heatmaps=args.boxAroundHeatmaps,
               dpi=args.dpi,
               interpolation_method=args.interpolationMethod)
