#!/usr/bin/env python
# -*- coding: utf-8 -*-


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
old_settings = np.seterr(all='ignore')
plt.ioff()


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        parents=[parserCommon.heatmapperMatrixArgs(),
                 parserCommon.heatmapperOutputArgs(mode='profile'),
                 parserCommon.heatmapperOptionalArgs(mode='profile')],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='This tool creates a profile plot for '
        'scores over sets of genomic regions. '
        'Typically, these regions are genes, but '
        'any other regions defined in BED '
        ' will work. A matrix generated '
        'by computeMatrix is required.',
        epilog='An example usage is: plotProfile -m <matrix file>',
        add_help=False)

    return parser


def process_args(args=None):
    args = parse_arguments().parse_args(args)

    # Because of galaxy, the value of this variables is normally
    # set to ''. Therefore this check is needed
    for attr in ['yMax', 'yMin']:
        try:
            args.__setattr__(attr, float(args.__getattribute__(attr)))
        except:
            args.__setattr__(attr, None)

    if args.plotHeight < 0.5:
        args.plotHeight = 0.5
    elif args.plotHeight > 100:
        args.plotHeight = 100

    return args


class Profile(object):

    def __init__(self, hm, out_file_name,
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
                 plots_per_row=8,
                 dpi=200):
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
        self.hm = hm
        self.out_file_name = out_file_name
        self.plot_title = plot_title
        self.y_axis_label = y_axis_label
        self.y_min = y_min
        self.y_max = y_max
        self.averagetype = averagetype
        self.reference_point_label = reference_point_label
        self.start_label = start_label
        self.end_label = end_label
        self.plot_height = plot_height
        self.plot_width = plot_width
        self.per_group = per_group
        self.plot_type = plot_type
        self.image_format = image_format
        self.color_list = color_list
        self.legend_location = legend_location
        self.plots_per_row = plots_per_row
        self.dpi = dpi

        # decide how many plots are needed
        if self.per_group:
            self.numplots = self.hm.matrix.get_num_groups()
            self.numlines = self.hm.matrix.get_num_samples()
        else:
            self.numplots = self.hm.matrix.get_num_samples()
            self.numlines = self.hm.matrix.get_num_groups()

        if self.numplots > self.plots_per_row:
            rows = np.ceil(self.numplots / float(self.plots_per_row)).astype(int)
            cols = self.plots_per_row
        else:
            rows = 1
            cols = self.numplots
        self.grids = gridspec.GridSpec(rows, cols)

        plt.rcParams['font.size'] = 10.0
        self.font_p = FontProperties()
        self.font_p.set_size('small')

        # convert cm values to inches
        plot_height_inches = rows * self.cm2inch(self.plot_height)[0]
        self.fig = plt.figure(figsize=self.cm2inch(cols * self.plot_width, rows * self.plot_height))
        self.fig.suptitle(self.plot_title, y=(1 - (0.06 / plot_height_inches)))
        self.xticks, self.xtickslabel = getProfileTicks(self.hm, self.reference_point_label,
                                                        self.start_label, self.end_label)

    @staticmethod
    def cm2inch(*tupl):
        inch = 2.54
        if isinstance(tupl[0], tuple):
            return tuple(i / inch for i in tupl[0])
        else:
            return tuple(i / inch for i in tupl)

    def plot_hexbin(self):
        from matplotlib import cm
        cmap = cm.coolwarm
        cmap.set_bad('black')

        for plot in range(self.numplots):
            col = plot % self.plots_per_row
            row = int(plot / float(self.plots_per_row))

            # split the ax to make room for the colorbar and for each of the
            # groups
            sub_grid = gridspec.GridSpecFromSubplotSpec(self.numlines, 2, subplot_spec=self.grids[row, col],
                                                        width_ratios=[0.92, 0.08], wspace=0.05, hspace=0.1)

            ax = self.fig.add_subplot(sub_grid[0, 0])

            ax.tick_params(
                axis='y',
                which='both',
                left='off',
                right='off',
                labelleft='on')

            if self.per_group:
                title = self.hm.matrix.group_labels[plot]
            else:
                title = self.hm.matrix.sample_labels[plot]

            vmin = np.inf
            vmax = -np.inf
            for data_idx in range(self.numlines):
                # get the max and min
                if self.per_group:
                    _row, _col = plot, data_idx
                else:
                    _row, _col = data_idx, plot

                sub_matrix = self.hm.matrix.get_matrix(_row, _col)
                ma = sub_matrix['matrix']
                x_values = np.tile(np.arange(ma.shape[1]), (ma.shape[0], 1))
                img = ax.hexbin(x_values.flatten(), ma.flatten(), cmap=cmap, mincnt=1)
                _vmin, _vmax = img.get_clim()
                if _vmin < vmin:
                    vmin = _vmin
                if _vmax > vmax:
                    vmax = _vmax
            self.fig.delaxes(ax)

            # iterate again after having computed the vmin and vmax
            ax_list = []
            for data_idx in range(self.numlines)[::-1]:
                ax = self.fig.add_subplot(sub_grid[data_idx, 0])
                if data_idx == 0:
                    ax.set_title(title)
                if data_idx != self.numlines - 1:
                    plt.setp(ax.get_xticklabels(), visible=False)

                if self.per_group:
                    _row, _col = plot, data_idx
                else:
                    _row, _col = data_idx, plot

                sub_matrix = self.hm.matrix.get_matrix(_row, _col)

                if self.per_group:
                    label = sub_matrix['sample']
                else:
                    label = sub_matrix['group']

                ma = sub_matrix['matrix']
                ax.set_axis_bgcolor('black')
                x_values = np.tile(np.arange(ma.shape[1]), (ma.shape[0], 1))
                img = ax.hexbin(x_values.flatten(), ma.flatten(), cmap=cmap, mincnt=1, vmin=vmin, vmax=vmax)

                # remove the numbers of the y axis for all plots
                ax.axes.set_ylabel(label)

                ax_list.append(ax)

                lims = ax.get_ylim()
                if self.y_min is not None:
                    lims = (self.y_min, lims[1])
                if self.y_max is not None:
                    lims = (lims[0], self.y_max)
                if lims[0] >= lims[1]:
                    lims = (lims[0], lims[0] + 1)
                ax.set_ylim(lims)

            ax_list[0].axes.set_xticks(self.xticks)
            ax_list[0].axes.set_xticklabels(self.xtickslabel)
            # align the first and last label
            # such that they don't fall off
            # the heatmap sides
            ticks = ax_list[-1].xaxis.get_major_ticks()
            ticks[0].label1.set_horizontalalignment('left')
            ticks[-1].label1.set_horizontalalignment('right')

            cax = self.fig.add_subplot(sub_grid[:, 1])
            self.fig.colorbar(img, cax=cax)

        plt.subplots_adjust(wspace=0.05, hspace=0.3)
        plt.tight_layout()
        plt.savefig(self.out_file_name, dpi=self.dpi, format=self.image_format)
        plt.close()

    def plot_heatmap(self):
        matrix_flatten = None
        if self.y_min is None:
            matrix_flatten = self.hm.matrix.flatten()
            # try to avoid outliers by using np.percentile
            self.y_min = np.percentile(matrix_flatten, 1.0)
            if np.isnan(self.y_min):
                self.y_min = None

        if self.y_max is None:
            if matrix_flatten is None:
                matrix_flatten = self.hm.matrix.flatten()
            # try to avoid outliers by using np.percentile
            self.y_max = np.percentile(matrix_flatten, 98.0)
            if np.isnan(self.y_max):
                self.y_max = None

        ax_list = []
        # turn off y ticks

        for plot in range(self.numplots):
            labels = []
            col = plot % self.plots_per_row
            row = int(plot / float(self.plots_per_row))

            # split the ax to make room for the colorbar
            sub_grid = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=self.grids[row, col],
                                                        width_ratios=[0.92, 0.08], wspace=0.05)

            ax = self.fig.add_subplot(sub_grid[0])
            cax = self.fig.add_subplot(sub_grid[1])

            ax.tick_params(
                axis='y',
                which='both',
                left='off',
                right='off',
                labelleft='on')

            if self.per_group:
                title = self.hm.matrix.group_labels[plot]
            else:
                title = self.hm.matrix.sample_labels[plot]

            ax.set_title(title)
            mat = []  # when drawing a heatmap (in contrast to drawing lines)
            for data_idx in range(self.numlines):
                if self.per_group:
                    row, col = plot, data_idx
                else:
                    row, col = data_idx, plot

                sub_matrix = self.hm.matrix.get_matrix(row, col)

                if self.per_group:
                    label = sub_matrix['sample']
                else:
                    label = sub_matrix['group']
                labels.append(label)

                mat.append(np.__getattribute__(self.averagetype)(sub_matrix['matrix'], axis=0))

            img = ax.imshow(np.vstack(mat), interpolation='nearest',
                            cmap='RdYlBu_r', aspect='auto', vmin=self.y_min, vmax=self.y_max)
            self.fig.colorbar(img, cax=cax)

            ax.axes.set_xticks(self.xticks)
            ax.axes.set_xticklabels(self.xtickslabel)
            # align the first and last label
            # such that they don't fall off
            # the heatmap sides
            ticks = ax.xaxis.get_major_ticks()
            ticks[0].label1.set_horizontalalignment('left')
            ticks[-1].label1.set_horizontalalignment('right')

            # add labels as y ticks labels
            ymin, ymax = ax.axes.get_ylim()
            pos, distance = np.linspace(ymin, ymax, len(labels), retstep=True, endpoint=False)
            d_half = float(distance) / 2
            yticks = [x + d_half for x in pos]

            ax.axes.set_yticks(yticks)
            # TODO: make rotation a parameter
            # ax.axes.set_yticklabels(labels[::-1], rotation='vertical')
            ax.axes.set_yticklabels(labels[::-1])

            ax_list.append(ax)

        plt.subplots_adjust(wspace=0.05, hspace=0.3)
        plt.tight_layout()
        plt.savefig(self.out_file_name, dpi=self.dpi, format=self.image_format)
        plt.close()

    def plot_profile(self):

        if not self.color_list:
            cmap_plot = plt.get_cmap('jet')
            if self.numlines > 1:
                # kmeans, so we need to color by cluster
                self.color_list = cmap_plot(np.arange(self.numlines, dtype=float) / float(self.numlines))
            else:
                self.color_list = cmap_plot(np.arange(self.numplots, dtype=float) / float(self.numplots))
        if (self.numlines > 1 and len(self.color_list) < self.numlines) or\
           (self.numlines == 1 and len(self.color_list) < self.numplots):
            sys.exit("\nThe given list of colors is too small, "
                     "at least {} colors are needed\n".format(self.numlines))
        for color in self.color_list:
            if not pltcolors.is_color_like(color):
                sys.exit("\nThe color name {} is not valid. Check "
                         "the name or try with a html hex string "
                         "for example #eeff22".format(color))

        first = True
        ax_list = []
        for plot in range(self.numplots):
            col = plot % self.plots_per_row
            row = int(plot / float(self.plots_per_row))
            if row == 0 and col == 0:
                ax = self.fig.add_subplot(self.grids[row, col])
            else:
                ax = self.fig.add_subplot(self.grids[row, col], sharey=ax_list[0])

            if self.per_group:
                title = self.hm.matrix.group_labels[plot]
                if row != 0:
                    plt.setp(ax.get_yticklabels(), visible=False)
            else:
                title = self.hm.matrix.sample_labels[plot]
                if col != 0:
                    plt.setp(ax.get_yticklabels(), visible=False)

            ax.set_title(title)
            for data_idx in range(self.numlines):
                if self.per_group:
                    _row, _col = plot, data_idx
                else:
                    _row, _col = data_idx, plot

                sub_matrix = self.hm.matrix.get_matrix(_row, _col)

                if self.per_group:
                    label = sub_matrix['sample']
                else:
                    label = sub_matrix['group']

                if self.numlines > 1:
                    coloridx = data_idx
                else:
                    coloridx = plot
                plot_single(ax, sub_matrix['matrix'],
                            self.averagetype,
                            self.color_list[coloridx],
                            label,
                            plot_type=self.plot_type)

            # remove the numbers of the y axis for all plots
            plt.setp(ax.get_yticklabels(), visible=False)

            if col == 0:
                # add the y axis label for the first plot
                # on each row and make the numbers and ticks visible
                plt.setp(ax.get_yticklabels(), visible=True)
                ax.axes.set_ylabel(self.y_axis_label)
                """
                # reduce the number of yticks by half
                num_ticks = len(ax.get_yticks())
                yticks = [ax.get_yticks()[i] for i in range(1, num_ticks, 2)]
                ax.set_yticks(yticks)
                """

            ax.axes.set_xticks(self.xticks)
            ax.axes.set_xticklabels(self.xtickslabel)
            # align the first and last label
            # such that they don't fall off
            # the heatmap sides
            ticks = ax.xaxis.get_major_ticks()
            ticks[0].label1.set_horizontalalignment('left')
            ticks[-1].label1.set_horizontalalignment('right')

            if first and self.plot_type not in ['heatmap', 'overlapped_lines']:
                ax.legend(loc=self.legend_location.replace('-', ' '),
                          ncol=1, prop=self.font_p,
                          frameon=False, markerscale=0.5)
                first = False

            """
            ax.legend(bbox_to_anchor=(-0.05, -1.13, 1., 1),
                      loc='upper center',
                      ncol=1, mode="expand", prop=font_p,
                      frameon=False, markerscale=0.5)
            """

            lims = ax.get_ylim()
            if self.y_min is not None:
                lims = (self.y_min, lims[1])
            if self.y_max is not None:
                lims = (lims[0], self.y_max)
            if lims[0] >= lims[1]:
                lims = (lims[0], lims[0] + 1)
            ax.set_ylim(lims)

            ax_list.append(ax)

        plt.subplots_adjust(wspace=0.05, hspace=0.3)
        plt.tight_layout()
        plt.savefig(self.out_file_name, dpi=self.dpi, format=self.image_format)
        plt.close()


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

    group_len_ratio = np.diff(hm.matrix.group_boundaries) / float(len(hm.matrix.regions))
    if np.any(group_len_ratio < 5.0 / 1000):
        problem = np.flatnonzero(group_len_ratio < 5.0 / 1000)
        sys.stderr.write("WARNING: Group '{}' is too small for plotting, you might want to remove it. \n".format(hm.matrix.group_labels[problem[0]]))

    if args.regionsLabel:
        hm.matrix.set_group_labels(args.regionsLabel)

    if args.samplesLabel and len(args.samplesLabel):
        hm.matrix.set_sample_labels(args.samplesLabel)

    if args.outFileNameData:
        hm.save_tabulated_values(args.outFileNameData, reference_point_label=args.refPointLabel,
                                 start_label=args.startLabel,
                                 end_label=args.endLabel,
                                 averagetype=args.averageType)

    if args.outFileSortedRegions:
        hm.save_BED(args.outFileSortedRegions)

    prof = Profile(hm, args.outFileName,
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
                   plots_per_row=args.numPlotsPerRow,
                   dpi=args.dpi)

    if args.plotType == 'heatmap':
        prof.plot_heatmap()
    elif args.plotType == 'overlapped_lines':
        prof.plot_hexbin()
    else:
        prof.plot_profile()
