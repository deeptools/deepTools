#!/usr/bin/env python
#-*- coding: utf-8 -*-

from __future__ import division
import sys

import argparse
import numpy as np
from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib import colors as pltcolors
import matplotlib.gridspec as gridspec

# own modules
from deeptools import parserCommon
from deeptools import heatmapper
from deeptools.heatmapper_utilities import plot_single, getProfileTicks

debug = 0
plt.ioff()


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        parents=[parserCommon.heatmapperMatrixArgs(),
                 parserCommon.heatmapperOutputArgs(mode='profile'),
                 parserCommon.heatmapperOptionalArgs(mode='profile')],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='This tool creates a profile plot for a '
        'score associated to genomic regions. '
        'Typically, these regions are genes, but '
        'any other regions defined in a BED or GFF '
        'format will work. A preprocessed matrix generated '
        'by the tool computeMatrix is required.',
        epilog='An example usage is: %(prog)s -m <matrix file>',
        add_help=False)

    return parser


def process_args(args=None):
    args = parse_arguments().parse_args(args)

    # Because of galaxy, the value of this variables is normally
    # set to ''. Therefore this check is needed
    for attr in ['yMax', 'yMin']:
        try:
            args.__setattr__(attr, float(args.__getattribute__(attr)))
#       except ValueError, TypeError:
        except:
            args.__setattr__(attr, None)

    if args.plotHeight < 0.5:
        args.plotHeight = 0.5
    elif args.plotHeight > 100:
        args.plotHeight = 100

    if args.regionsLabel != 'genes':
        args.regionsLabel = \
            [x.strip() for x in args.regionsLabel.split(',')]

        if len(set(args.regionsLabel)) != len(args.regionsLabel):
            print "The group labels given contain repeated names. Please "
            "give a unique name to each value. The values given are "
            "{}\n".format(args.regionsLabel)
            exit(1)
    else:
        args.regionsLabel = []

    return args


def plot_profile(hm, out_file_name,
                 plot_title='', y_axis_label='',
                 y_min=None, y_max=None,
                 averagetype='median',
                 reference_point_label='TSS',
                 start_label='TSS', end_label="TES",
                 plot_height=7,
                 plot_width=11,
                 per_group=False,
                 plot_type='simple',
                 image_format=None,
                 color_list=None,
                 legend_location='auto',
                 plots_per_row=8):
    """
    Using the hm matrix, makes a line plot
    either per group or per sample
    using the specified parameters.

    Args:
        hm: heatmapper object
        out_file_name: string
        plot_title: string
        y_axis_label: list
        y_min: int
        y_max: int
        averagetype: mean, sum, median
        reference_point_label: string
        start_label: string
        end_label: string
        plot_height: in cm
        plot_width: in cm
        per_group: bool
        plot_type: string
        image_format: string
        color_list: list
        legend_location:
        plots_per_row: int

    Returns:

    """
    # the following line is to temporary compute the log
    # of the matrix
    # hm.matrix.matrix = np.log(hm.matrix.matrix)

    # decide how many plots are needed
    if per_group:
        numplots = hm.matrix.get_num_groups()
        numlines = hm.matrix.get_num_samples()
    else:
        numplots = hm.matrix.get_num_samples()
        numlines = hm.matrix.get_num_groups()

    if numplots > plots_per_row:
        rows = np.ceil(numplots / plots_per_row).astype(int)
        cols = plots_per_row
    else:
        rows = 1
        cols = numplots
    grids = gridspec.GridSpec(rows, cols)

    plt.rcParams['font.size'] = 10.0
    font_p = FontProperties()
    font_p.set_size('small')
#    rcParams['font.size'] = 9.0

    # convert cm values to inches
    plot_height_inches = rows * float(plot_height) / 2.54
    fig_width = cols * float(plot_width) / 2.54
    fig = plt.figure(figsize=(fig_width, plot_height_inches))

    # add xticks and labels to the last plot
    # define the xticks

    fig.suptitle(plot_title, y=(1 - (0.06 / plot_height_inches)))

    if not color_list:
        cmap_plot = plt.get_cmap('jet')
        if numlines > 1:
            #kmeans, so we need to color by cluster
            color_list = cmap_plot(np.arange(numlines, dtype=float) / numlines)
        else :
            color_list = cmap_plot(np.arange(numplots, dtype=float) / numplots)
    else:
        if (numlines > 1 and len(color_list) < numlines) or (numlines == 1 and len(color_list) < numplots) :
            sys.stderr.write("\nThe given list of colors is too small, "
                             "at least {} colors are needed\n".format(numlines))
            exit(1)
        for color in color_list:
            if not pltcolors.is_color_like(color):
                sys.stderr.write("\nThe color name {} is not valid. Check "
                                 "the name or try with a html hex string "
                                 "for example #eeff22".format(color))

                exit(1)

    #If we plot a heatmap, we need to normalize the values to be between 0 and 1
    if plot_type == 'heatmap':
        max_val = None
        for plot in range(numplots):
            col = plot % plots_per_row
            row = int(plot / plots_per_row)
            for data_idx in range(numlines):
                if per_group:
                    row, col = plot, data_idx
                else:
                    row, col = data_idx, plot
                sub_matrix = hm.matrix.get_matrix(row, col)
                cur_max = np.amax(np.__getattribute__(averagetype)(sub_matrix['matrix'], axis=0))
                if max_val is None or cur_max > max_val:
                    max_val = cur_max

    xticks, xtickslabel = getProfileTicks(hm, reference_point_label,
                                          start_label, end_label)
    first = True
    sample_ymax = None
    sample_ymin = None
    ax_list = []
    for plot in range(numplots):
        col = plot % plots_per_row
        row = int(plot / plots_per_row)
        if row == 0 and col == 0:
            ax = fig.add_subplot(grids[row, col])
        else:
            ax = fig.add_subplot(grids[row, col], sharey=ax_list[0])

        if per_group:
            title = hm.matrix.group_labels[plot]
            if row != 0:
                plt.setp(ax.get_yticklabels(), visible=False)
        else:
            title = hm.matrix.sample_labels[plot]
            if col != 0:
                plt.setp(ax.get_yticklabels(), visible=False)

        ax.set_title(title)
        mat = []  # when drawing a heatmap (in contrast to drawing lines)
        for data_idx in range(numlines):
            if per_group:
                row, col = plot, data_idx
            else:
                row, col = data_idx, plot

            sub_matrix = hm.matrix.get_matrix(row, col)

            if per_group:
                label = sub_matrix['sample']
            else:
                label = sub_matrix['group']

            if plot_type == 'heatmap':
                # if plotting a heatmap, the individual rows
                # need to be collected before plotting
                mat.append(np.__getattribute__(averagetype)(sub_matrix['matrix'],
                                                            axis=0))
            else:
                if numlines > 1:
                    coloridx = data_idx
                else :
                    coloridx = plot
                plot_single(ax, sub_matrix['matrix'],
                            averagetype,
                            color_list[coloridx],
                            label,
                            plot_type=plot_type)

        if plot_type == 'heatmap':
            ax.imshow(np.divide(np.asmatrix(mat), max_val), interpolation='nearest', cmap='seismic_r')

        if (per_group and row > 0) or (per_group is False and col > 0):
            # remove the numbers of the y axis for all plots
            # except the first one
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            # add the y axis label for the first plot
            ax.axes.set_ylabel(y_axis_label)
            """
            # reduce the number of yticks by half
            num_ticks = len(ax.get_yticks())
            yticks = [ax.get_yticks()[i] for i in range(1, num_ticks, 2)]
            ax.set_yticks(yticks)
            """

        ax.axes.set_xticks(xticks)
        ax.axes.set_xticklabels(xtickslabel)
        # align the first and last label
        # such that they don't fall off 
        # the heatmap sides
        ticks = ax.xaxis.get_major_ticks()
        ticks[0].label1.set_horizontalalignment('left')
        ticks[-1].label1.set_horizontalalignment('right')

        if first:
            ax.legend(loc=legend_location.replace('-', ' '),
                      ncol=1, prop=font_p,
                      frameon=False, markerscale=0.5)
            first = False
        """
        ax.legend(bbox_to_anchor=(-0.05, -1.13, 1., 1),
                  loc='upper center',
                  ncol=1, mode="expand", prop=font_p,
                  frameon=False, markerscale=0.5)
        """

        if plot_type != 'heatmap':
            lims = ax.get_ylim()
            if y_min is not None:
                lims = (y_min, lims[1])
            if y_max is not None:
                lims = (lims[0], y_max)
            if lims[0] >= lims[1]:
                lims = (lims[0], lims[0]+1)
            ax.set_ylim(lims)

        ax_list.append(ax)

    plt.subplots_adjust(wspace=0.05, hspace=0.3)
    plt.tight_layout()
    plt.savefig(out_file_name, dpi=200, format=image_format)


def main(args=None):
    args = process_args(args)
    hm = heatmapper.heatmapper()
    matrix_file = args.matrixFile.name
    args.matrixFile.close()
    hm.read_matrix_file(matrix_file, default_group_name=args.regionsLabel)

    if args.kmeans is not None:
        hm.matrix.hmcluster(args.kmeans, method='kmeans')

    if len(args.regionsLabel):
        hm.matrix.set_group_labels(args.regionsLabel)

    if args.samplesLabel and len(args.samplesLabel):
        hm.matrix.set_sample_labels(args.samplesLabel)

    #if args.outFileNameData:
    #    hm.saveTabulatedValues(args.outFileNameData)

    if args.outFileSortedRegions:
        hm.save_BED(args.outFileSortedRegions)

    plot_profile(hm, args.outFileName,
                 plot_title=args.plotTitle,
                 y_axis_label=args.yAxisLabel,
                 y_min=args.yMin, y_max=args.yMax,
                 averagetype=args.averageType,
                 reference_point_label=args.refPointLabel,
                 start_label=args.startLabel,
                 end_label=args.endLabel,
                 plot_height=args.plotHeight,
                 plot_width=args.plotWidth,
                 per_group=args.perGroup,
                 plot_type=args.plotType,
                 image_format=args.plotFileFormat,
                 color_list=args.colors,
                 legend_location=args.legendLocation,
                 plots_per_row=args.numPlotsPerRow)
