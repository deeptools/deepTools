#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys

import argparse
import numpy as np
from math import ceil
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
import deeptools.cm  # noqa: F401
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib import colors as pltcolors
import matplotlib.gridspec as gridspec

import plotly.offline as py
import plotly.graph_objs as go

# own modules
from deeptools import parserCommon
from deeptools import heatmapper
from deeptools.heatmapper_utilities import plot_single, plotly_single, getProfileTicks
from deeptools.computeMatrixOperations import filterHeatmapValues


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

    # Ensure that yMin/yMax are there and a list
    try:
        assert args.yMin is not None
    except:
        args.yMin = [None]
    try:
        assert args.yMax is not None
    except:
        args.yMax = [None]

    # Sometimes Galaxy sends --yMax '' and --yMin ''
    if args.yMin == ['']:
        args.yMin = [None]
    if args.yMax == ['']:
        args.yMax = [None]

    # Convert to floats
    if args.yMin != [None]:
        foo = [float(x) for x in args.yMin]
        args.yMin = foo
    if args.yMax != [None]:
        foo = [float(x) for x in args.yMax]
        args.yMax = foo

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
                 reference_point_label=None,
                 start_label='TSS', end_label='TES',
                 plot_height=7,
                 plot_width=11,
                 per_group=False,
                 plot_type='lines',
                 image_format=None,
                 color_list=None,
                 legend_location='best',
                 plots_per_row=8,
                 label_rotation=0,
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
            y_min: list
            y_max: list
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
            label_rotation: float

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
        self.label_rotation = label_rotation
        self.dpi = dpi

        # Honor reference point labels from computeMatrix
        if reference_point_label is None:
            self.reference_point_label = hm.parameters['ref point']

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

        # Ensure that the labels are vectors
        nSamples = len(self.hm.matrix.sample_labels)
        if not isinstance(self.reference_point_label, list):
            self.reference_point_label = [self.reference_point_label] * nSamples
        if not isinstance(self.start_label, list):
            self.start_label = [self.start_label] * nSamples
        if not isinstance(self.end_label, list):
            self.end_label = [self.end_label] * nSamples

    def getTicks(self, idx):
        """
        This is essentially a wrapper around getProfileTicks to accomdate the fact that each column has its own ticks.
        """
        xticks, xtickslabel = getProfileTicks(self.hm, self.reference_point_label[idx], self.start_label[idx], self.end_label[idx], idx)
        return xticks, xtickslabel

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

        if self.image_format == "plotly":
            return self.plotly_hexbin()

        for plot in range(self.numplots):
            col = plot % self.plots_per_row
            row = int(plot / float(self.plots_per_row))
            localYMin = None
            localYMax = None

            # split the ax to make room for the colorbar and for each of the
            # groups
            sub_grid = gridspec.GridSpecFromSubplotSpec(self.numlines, 2, subplot_spec=self.grids[row, col],
                                                        width_ratios=[0.92, 0.08], wspace=0.05, hspace=0.1)

            ax = self.fig.add_subplot(sub_grid[0, 0])

            ax.tick_params(
                axis='y',
                which='both',
                left=False,
                right=False,
                labelleft=True)

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

                if localYMin is None or self.y_min[col % len(self.y_min)] < localYMin:
                    localYMin = self.y_min[col % len(self.y_min)]
                if localYMax is None or self.y_max[col % len(self.y_max)] > localYMax:
                    localYMax = self.y_max[col % len(self.y_max)]
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
                try:
                    # matplotlib 2.0
                    ax.set_facecolor('black')
                except:
                    # matplotlib <2.0
                    ax.set_axis_bgcolor('black')
                x_values = np.tile(np.arange(ma.shape[1]), (ma.shape[0], 1))
                img = ax.hexbin(x_values.flatten(), ma.flatten(), cmap=cmap, mincnt=1, vmin=vmin, vmax=vmax)

                if plot == 0:
                    ax.axes.set_ylabel(label)

                ax_list.append(ax)

                lims = ax.get_ylim()
                if localYMin is not None:
                    lims = (localYMin, lims[1])
                if localYMax is not None:
                    lims = (lims[0], localYMax)
                if lims[0] >= lims[1]:
                    lims = (lims[0], lims[0] + 1)
                ax.set_ylim(lims)

            xticks, xtickslabel = self.getTicks(plot)
            if np.ceil(max(xticks)) != float(ma.shape[1] - 1):
                tickscale = float(sub_matrix['matrix'].shape[1]) / max(xticks)
                xticks_use = [x * tickscale for x in xticks]
                ax_list[0].axes.set_xticks(xticks_use)
            else:
                ax_list[0].axes.set_xticks(xticks)
            ax_list[0].axes.set_xticklabels(xtickslabel, rotation=self.label_rotation)
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

    def plotly_hexbin(self):
        """plot_hexbin, but for plotly. it's annoying that we have to have sub-subplots"""
        fig = go.Figure()
        cols = self.plots_per_row if self.numplots > self.plots_per_row else self.numplots
        rows = np.ceil(self.numplots / float(cols)).astype(int)
        fig['layout'].update(title=self.plot_title)
        domainWidth = .9 / cols
        domainHeight = .9 / rows
        bufferHeight = 0.0
        if rows > 1:
            bufferHeight = 0.1 / (rows - 1)
        else:
            domainHeight = 1.0
        bufferWidth = 0.0
        if cols > 1:
            bufferWidth = 0.1 / (cols - 1)
        else:
            domainWidth = 1.0
        subHeight = domainHeight / float(self.numlines)
        if self.per_group:
            sideLabels = self.hm.matrix.sample_labels
        else:
            sideLabels = self.hm.matrix.group_labels

        data = []
        annos = []
        vmin = np.inf
        vmax = -np.inf
        for i in range(self.numplots):
            row = rows - i / self.plots_per_row - 1
            col = i % self.plots_per_row

            if self.per_group:
                title = self.hm.matrix.group_labels[i]
            else:
                title = self.hm.matrix.sample_labels[i]

            base = row * (domainHeight + bufferHeight)
            domain = [base, base + domainHeight]
            titleY = base + domainHeight
            base = col * (domainWidth + bufferWidth)
            domain = [base, base + domainWidth]
            titleX = base + 0.5 * domainWidth
            xanchor = 'x{}'.format(i + 1)
            fig['layout']['xaxis{}'.format(i + 1)] = dict(domain=domain)
            annos.append({'yanchor': 'bottom', 'xref': 'paper', 'xanchor': 'center', 'yref': 'paper', 'text': title, 'y': titleY, 'x': titleX, 'font': {'size': 16}, 'showarrow': False})

            # set yMin/yMax
            yMin = np.inf
            yMax = -np.inf
            for j in range(self.numlines):
                # get the max and min
                if self.per_group:
                    _row, _col = i, j
                else:
                    _row, _col = j, i

                ma = self.hm.matrix.get_matrix(_row, _col)['matrix']
                if np.min(ma) < yMin:
                    yMin = np.min(ma)
                if np.max(ma) > yMax:
                    yMax = np.max(ma)
            if self.y_min[i % len(self.y_min)] is not None:
                yMin = self.y_min[i % len(self.y_min)]
            if self.y_max[i % len(self.y_max)] is not None:
                yMax = self.y_max[i % len(self.y_max)]

            for j in range(self.numlines):
                if self.per_group:
                    _row, _col = i, j
                else:
                    _row, _col = j, i
                foo = i * self.numlines + j + 1
                yanchor = 'y{}'.format(foo)
                base = row * (domainHeight + bufferHeight) + j * subHeight
                domain = [base, base + subHeight]
                fig['layout']['yaxis{}'.format(foo)] = {'domain': domain, 'title': self.y_axis_label, 'anchor': xanchor, 'range': [yMin, yMax]}
                if j == 0:
                    _ = "xaxis{}".format(xanchor[1:])
                    fig['layout'][_].update(anchor='y{}'.format(foo))
                if col == 0:
                    titleY = base + 0.5 * subHeight
                    annos.append({'yanchor': 'middle', 'xref': 'paper', 'xanchor': 'left', 'yref': 'paper', 'text': sideLabels[j], 'y': titleY, 'x': -0.03, 'font': {'size': 16}, 'showarrow': False, 'textangle': -90})

                sub_matrix = self.hm.matrix.get_matrix(_row, _col)
                ma = self.hm.matrix.get_matrix(_row, _col)['matrix']

                fig['layout']['xaxis{}'.format(i + 1)].update(range=[0, ma.shape[1]])

                if self.per_group:
                    label = sub_matrix['sample']
                else:
                    label = sub_matrix['group']

                # Manually compute the 2D histogram with 100x100 bins
                x_values = np.tile(np.arange(ma.shape[1]), (ma.shape[0], 1))
                z, xe, ye = np.histogram2d(x_values.flatten(), ma.flatten(), bins=100, range=[[0, ma.shape[1]], [yMin, yMax]])

                _vmin = np.min(z)
                _vmax = np.max(z)
                if _vmin < vmin:
                    vmin = _vmin
                if _vmax > vmax:
                    vmax = _vmax

                trace = go.Contour(z=z.T, x=xe, y=ye, xaxis=xanchor, yaxis=yanchor, name=label, connectgaps=False)
                data.append(trace)

            # Assume the bounds for the last graph are correct
            totalWidth = ma.shape[1]
            xticks, xtickslabel = self.getTicks(i)
            if np.ceil(max(xticks)) != float(totalWidth):
                tickscale = float(totalWidth) / max(xticks)
                xticks_use = [x * tickscale for x in xticks]
            else:
                xticks_use = xticks
            xticks_use = [np.ceil(x) for x in xticks_use]
            fig['layout']['xaxis{}'.format(i + 1)].update(tickmode='array', tickvals=xticks_use, ticktext=xtickslabel, tickangle=self.label_rotation)

        for trace in data:
            trace.update(zmin=vmin, zmax=vmax)

        fig.add_traces(data)
        fig['layout']['annotations'] = annos
        py.plot(fig, filename=self.out_file_name, auto_open=False)

    def plot_heatmap(self):
        cmap = ['RdYlBu_r']
        if self.color_list is not None:  # check the length to be equal to the numebr of plots otherwise multiply it!
            cmap = self.color_list
        if len(cmap) < self.numplots:
            all_colors = cmap
            for i in range(ceil(self.numplots / len(cmap))):
                cmap.extend(all_colors)
        matrix_flatten = None
        if self.y_min == [None]:
            matrix_flatten = self.hm.matrix.flatten()
            # try to avoid outliers by using np.percentile
            self.y_min = [np.percentile(matrix_flatten, 1.0)]
            if np.isnan(self.y_min[0]):
                self.y_min = [None]

        if self.y_max == [None]:
            if matrix_flatten is None:
                matrix_flatten = self.hm.matrix.flatten()
            # try to avoid outliers by using np.percentile
            self.y_max = [np.percentile(matrix_flatten, 98.0)]
            if np.isnan(self.y_max[0]):
                self.y_max = [None]

        if self.image_format == "plotly":
            return self.plotly_heatmap()

        ax_list = []
        # turn off y ticks
        for plot in range(self.numplots):
            labels = []
            col = plot % self.plots_per_row
            row = int(plot / float(self.plots_per_row))
            localYMin = None
            localYMax = None

            # split the ax to make room for the colorbar
            sub_grid = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=self.grids[row, col],
                                                        width_ratios=[0.92, 0.08], wspace=0.05)

            ax = self.fig.add_subplot(sub_grid[0])
            cax = self.fig.add_subplot(sub_grid[1])

            ax.tick_params(
                axis='y',
                which='both',
                left=False,
                right=False,
                labelleft=True)

            if self.per_group:
                title = self.hm.matrix.group_labels[plot]
                tickIdx = plot % self.hm.matrix.get_num_samples()
            else:
                title = self.hm.matrix.sample_labels[plot]
                tickIdx = plot
            ax.set_title(title)
            mat = []  # when drawing a heatmap (in contrast to drawing lines)
            for data_idx in range(self.numlines):
                if self.per_group:
                    row, col = plot, data_idx
                else:
                    row, col = data_idx, plot
                if localYMin is None or self.y_min[col % len(self.y_min)] < localYMin:
                    localYMin = self.y_min[col % len(self.y_min)]
                if localYMax is None or self.y_max[col % len(self.y_max)] > localYMax:
                    localYMax = self.y_max[col % len(self.y_max)]

                sub_matrix = self.hm.matrix.get_matrix(row, col)

                if self.per_group:
                    label = sub_matrix['sample']
                else:
                    label = sub_matrix['group']
                labels.append(label)
                mat.append(np.ma.__getattribute__(self.averagetype)(sub_matrix['matrix'], axis=0))
            img = ax.imshow(np.vstack(mat), interpolation='nearest',
                            cmap=cmap[plot], aspect='auto', vmin=localYMin, vmax=localYMax)
            self.fig.colorbar(img, cax=cax)

            totalWidth = np.vstack(mat).shape[1]
            xticks, xtickslabel = self.getTicks(tickIdx)
            if np.ceil(max(xticks)) != float(totalWidth - 1):
                tickscale = float(totalWidth) / max(xticks)
                xticks_use = [x * tickscale for x in xticks]
                ax.axes.set_xticks(xticks_use)
            else:
                ax.axes.set_xticks(xticks)
            ax.axes.set_xticklabels(xtickslabel, rotation=self.label_rotation)
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

            # TODO: make rotation a parameter
            # ax.axes.set_yticklabels(labels[::-1], rotation='vertical')
            if plot == 0:
                ax.axes.set_yticks(yticks)
                ax.axes.set_yticklabels(labels[::-1])
            else:
                ax.axes.set_yticklabels([])
            # matplotlib 3.1.1 (and likely some earlier versions) will change the ylim if you change the tick locations!
            ax.axes.set_ylim([ymin, ymax])

            ax_list.append(ax)

        plt.subplots_adjust(wspace=0.05, hspace=0.3)
        plt.tight_layout()
        plt.savefig(self.out_file_name, dpi=self.dpi, format=self.image_format)
        plt.close()

    def plotly_heatmap(self):
        """plot_heatmap, but with plotly output"""
        fig = go.Figure()
        cols = self.plots_per_row if self.numplots > self.plots_per_row else self.numplots
        rows = np.ceil(self.numplots / float(cols)).astype(int)
        fig['layout'].update(title=self.plot_title)
        domainWidth = .9 / cols
        domainHeight = .9 / rows
        bufferHeight = 0.0
        if rows > 1:
            bufferHeight = 0.1 / (rows - 1)
        else:
            domainHeight = 1.0
        bufferWidth = 0.0
        if cols > 1:
            bufferWidth = 0.1 / (cols - 1)
        else:
            domainWidth = 1.0

        data = []
        annos = []
        zmin = np.inf
        zmax = -np.inf
        for i in range(self.numplots):
            row = rows - i / self.plots_per_row - 1
            col = i % self.plots_per_row

            if self.per_group:
                title = self.hm.matrix.group_labels[i]
            else:
                title = self.hm.matrix.sample_labels[i]

            base = row * (domainHeight + bufferHeight)
            domain = [base, base + domainHeight]
            titleY = base + domainHeight
            xanchor = 'x{}'.format(i + 1)
            yanchor = 'y{}'.format(i + 1)
            visible = False
            if col == 0:
                visible = True
            fig['layout']['yaxis{}'.format(i + 1)] = {'domain': domain, 'anchor': xanchor, 'visible': visible}
            base = col * (domainWidth + bufferWidth)
            domain = [base, base + domainWidth]
            titleX = base + 0.5 * domainWidth
            fig['layout']['xaxis{}'.format(i + 1)] = {'domain': domain, 'anchor': yanchor}
            annos.append({'yanchor': 'bottom', 'xref': 'paper', 'xanchor': 'center', 'yref': 'paper', 'text': title, 'y': titleY, 'x': titleX, 'font': {'size': 16}, 'showarrow': False})

            mat = []
            labels = []
            for j in range(self.numlines):
                if self.per_group:
                    row, col = i, j
                else:
                    row, col = j, i

                sub_matrix = self.hm.matrix.get_matrix(row, col)

                if self.per_group:
                    label = sub_matrix['sample']
                else:
                    label = sub_matrix['group']
                labels.append(label)
                mat.append(np.ma.__getattribute__(self.averagetype)(sub_matrix['matrix'], axis=0))
                if np.min(mat[-1]) < zmin:
                    zmin = np.min(mat[-1])
                if np.max(mat[-1]) > zmax:
                    zmax = np.max(mat[-1])
            totalWidth = len(mat[-1])
            trace = go.Heatmap(name=title, z=mat, x=range(totalWidth + 1), y=labels, xaxis=xanchor, yaxis=yanchor)
            data.append(trace)

            # Add ticks
            xticks, xtickslabel = self.getTicks(i)
            if np.ceil(max(xticks)) != float(totalWidth):
                tickscale = float(totalWidth) / max(xticks)
                xticks_use = [x * tickscale for x in xticks]
            else:
                xticks_use = xticks
            xticks_use = [np.ceil(x) for x in xticks_use]
            fig['layout']['xaxis{}'.format(i + 1)].update(tickmode='array', tickvals=xticks_use, ticktext=xtickslabel, tickangle=self.label_rotation)

        # Adjust color scale limits
        for i, trace in enumerate(data):
            zminUse = zmin
            zmaxUse = zmax
            if self.y_min[i % len(self.y_min)] is not None:
                zminUse = self.y_min[i % len(self.y_min)]
            if self.y_max[i % len(self.y_max)] is not None:
                zmaxUse = self.y_max[i % len(self.y_max)]
            trace.update(zmin=zminUse, zmax=zmaxUse)

        fig.add_traces(data)
        fig['layout']['annotations'] = annos
        py.plot(fig, filename=self.out_file_name, auto_open=False)

    def plot_profile(self):
        if self.y_min is None:
            self.y_min = [None]
        if self.y_max is None:
            self.y_max = [None]

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

        if self.image_format == "plotly":
            return self.plotly_profile()

        first = True
        ax_list = []
        globalYmin = np.inf
        globalYmax = -np.inf
        for plot in range(self.numplots):
            localYMin = None
            localYMax = None
            col = plot % self.plots_per_row
            row = int(plot / float(self.plots_per_row))
            if (row == 0 and col == 0) or len(self.y_min) > 1 or len(self.y_max) > 1:
                ax = self.fig.add_subplot(self.grids[row, col])
            else:
                ax = self.fig.add_subplot(self.grids[row, col])

            if self.per_group:
                title = self.hm.matrix.group_labels[plot]
                if row != 0 and len(self.y_min) == 1 and len(self.y_max) == 1:
                    plt.setp(ax.get_yticklabels(), visible=False)
                tickIdx = plot % self.hm.matrix.get_num_samples()
            else:
                title = self.hm.matrix.sample_labels[plot]
                if col != 0 and len(self.y_min) == 1 and len(self.y_max) == 1:
                    plt.setp(ax.get_yticklabels(), visible=False)
                tickIdx = plot

            ax.set_title(title)
            for data_idx in range(self.numlines):
                if self.per_group:
                    _row, _col = plot, data_idx
                else:
                    _row, _col = data_idx, plot
                if localYMin is None or self.y_min[col % len(self.y_min)] < localYMin:
                    localYMin = self.y_min[col % len(self.y_min)]
                if localYMax is None or self.y_max[col % len(self.y_max)] > localYMax:
                    localYMax = self.y_max[col % len(self.y_max)]

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
            globalYmin = min(np.float64(globalYmin), ax.get_ylim()[0])
            globalYmax = max(globalYmax, ax.get_ylim()[1])

            # Exclude ticks from all but one subplot by default
            if col > 0 and len(self.y_min) == 1 and len(self.y_max) == 1:
                plt.setp(ax.get_yticklabels(), visible=False)

            totalWidth = sub_matrix['matrix'].shape[1]
            xticks, xtickslabel = self.getTicks(tickIdx)
            if np.ceil(max(xticks)) != float(totalWidth - 1):
                tickscale = float(totalWidth) / max(xticks)
                xticks_use = [x * tickscale for x in xticks]
                ax.axes.set_xticks(xticks_use)
            else:
                ax.axes.set_xticks(xticks)
            ax.axes.set_xticklabels(xtickslabel, rotation=self.label_rotation)
            # align the first and last label
            # such that they don't fall off
            # the heatmap sides
            ticks = ax.xaxis.get_major_ticks()
            ticks[0].label1.set_horizontalalignment('left')
            ticks[-1].label1.set_horizontalalignment('right')

            if first and self.y_axis_label != '':
                ax.set_ylabel(self.y_axis_label)
            if first and self.plot_type not in ['heatmap', 'overlapped_lines']:
                ax.legend(loc=self.legend_location.replace('-', ' '),
                          ncol=1, prop=self.font_p,
                          frameon=False, markerscale=0.5)
                if len(self.y_min) == 1 and len(self.y_max) == 1:
                    first = False
            ax_list.append(ax)

        # It turns out that set_ylim only takes np.float64s
        for sample_id, subplot in enumerate(ax_list):
            localYMin = self.y_min[sample_id % len(self.y_min)]
            localYMax = self.y_max[sample_id % len(self.y_max)]
            lims = [globalYmin, globalYmax]
            if localYMin is not None:
                if localYMax is not None:
                    lims = (np.float64(localYMin), np.float64(localYMax))
                else:
                    lims = (np.float64(localYMin), lims[1])
            elif localYMax is not None:
                lims = (lims[0], np.float64(localYMax))
            if lims[0] >= lims[1]:
                lims = (lims[0], lims[0] + 1)
            ax_list[sample_id].set_ylim(lims)

        plt.subplots_adjust(wspace=0.05, hspace=0.3)
        plt.tight_layout()
        plt.savefig(self.out_file_name, dpi=self.dpi, format=self.image_format)
        plt.close()

    def plotly_profile(self):
        """
        plot_profile for plotly output

        y_min, y_max, and color_list are set already
        """
        fig = go.Figure()
        cols = self.plots_per_row if self.numplots > self.plots_per_row else self.numplots
        rows = np.ceil(self.numplots / float(cols)).astype(int)
        fig['layout'].update(title=self.plot_title)
        domainWidth = .9 / cols
        domainHeight = .9 / rows
        bufferHeight = 0.0
        if rows > 1:
            bufferHeight = 0.1 / (rows - 1)
        bufferWidth = 0.0
        if cols > 1:
            bufferWidth = 0.1 / (cols - 1)

        data = []
        annos = []
        yMin = None
        yMax = None
        for i in range(self.numplots):
            row = np.floor(i / self.plots_per_row)
            # row = rows - i / self.plots_per_row - 1
            col = i % self.plots_per_row
            xanchor = 'x{}'.format(i + 1)
            yanchor = 'y{}'.format(i + 1)
            base = row * (domainHeight + bufferHeight)
            domain = [base, base + domainHeight]
            titleY = base + domainHeight
            fig['layout']['yaxis{}'.format(i + 1)] = {'domain': domain, 'title': self.y_axis_label, 'anchor': xanchor, 'autorange': False}
            base = col * (domainWidth + bufferWidth)
            domain = [base, base + domainWidth]
            titleX = base + 0.5 * domainWidth
            fig['layout']['xaxis{}'.format(i + 1)] = {'domain': domain, 'anchor': yanchor}

            if self.per_group:
                title = self.hm.matrix.group_labels[i]
            else:
                title = self.hm.matrix.sample_labels[i]
            annos.append({'yanchor': 'bottom', 'xref': 'paper', 'xanchor': 'center', 'yref': 'paper', 'text': title, 'y': titleY, 'x': titleX, 'font': {'size': 16}, 'showarrow': False})

            for j in range(self.numlines):
                if self.per_group:
                    _row, _col = i, j
                else:
                    _row, _col = j, i

                sub_matrix = self.hm.matrix.get_matrix(_row, _col)
                fig['layout']['xaxis{}'.format(i + 1)].update(range=[0, sub_matrix['matrix'].shape[1]])

                if self.per_group:
                    label = sub_matrix['sample']
                else:
                    label = sub_matrix['group']

                if self.numlines > 1:
                    coloridx = j
                else:
                    coloridx = i
                color = self.color_list[coloridx]
                traces = plotly_single(sub_matrix['matrix'],
                                       self.averagetype,
                                       color,
                                       label,
                                       plot_type=self.plot_type)
                for trace in traces:
                    trace.update(xaxis=xanchor, yaxis=yanchor)
                    if yMin is None or min(trace['y']) < yMin:
                        yMin = min(trace['y'])
                    if yMax is None or max(trace['y']) > yMax:
                        yMax = max(trace['y'])
                if row == col == 0:
                    traces[0].update(showlegend=True)
                data.extend(traces)
            totalWidth = sub_matrix['matrix'].shape[1]
            xticks, xtickslabel = self.getTicks(i)
            if np.ceil(max(xticks)) != float(totalWidth):
                tickscale = float(totalWidth) / max(xticks)
                xticks_use = [x * tickscale for x in xticks]
            else:
                xticks_use = xticks
            xticks_use = [np.ceil(x) for x in xticks_use]
            fig['layout']['xaxis{}'.format(i + 1)].update(tickmode='array', tickvals=xticks_use, ticktext=xtickslabel, tickangle=self.label_rotation)

        # Set the y limits
        for i in range(self.numplots):
            yaxis = 'yaxis{}'.format(i + 1)
            yRange = [yMin, yMax]
            if self.y_min[i % len(self.y_min)] is not None:
                yRange[0] = self.y_min[i % len(self.y_min)]
            if self.y_max[i % len(self.y_max)] is not None:
                yRange[1] = self.y_max[i % len(self.y_max)]
            fig['layout'][yaxis].update(range=yRange)

        fig.add_traces(data)
        fig['layout']['annotations'] = annos
        py.plot(fig, filename=self.out_file_name, auto_open=False)


def main(args=None):
    args = process_args(args)
    hm = heatmapper.heatmapper()
    matrix_file = args.matrixFile.name
    args.matrixFile.close()
    hm.read_matrix_file(matrix_file)

    if hm.parameters['min threshold'] is not None or hm.parameters['max threshold'] is not None:
        filterHeatmapValues(hm, hm.parameters['min threshold'], hm.parameters['max threshold'])

    if args.kmeans is not None:
        hm.matrix.hmcluster(args.kmeans, method='kmeans', clustering_samples=args.clusterUsingSamples)
    else:
        if args.hclust is not None:
            print("Performing hierarchical clustering."
                  "Please note that it might be very slow for large datasets.\n")
            hm.matrix.hmcluster(args.hclust, method='hierarchical', clustering_samples=args.clusterUsingSamples)

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
                   label_rotation=args.label_rotation,
                   dpi=args.dpi)

    if args.plotType == 'heatmap':
        prof.plot_heatmap()
    elif args.plotType == 'overlapped_lines':
        prof.plot_hexbin()
    else:
        prof.plot_profile()
