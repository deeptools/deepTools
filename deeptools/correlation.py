import sys
import itertools
import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.stats
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
from deeptools import cm  # noqa: F401
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker
import matplotlib.mlab
import matplotlib.markers
import matplotlib.colors as pltcolors
from deeptools.utilities import toString, convertCmap

import plotly.offline as offline
import plotly.graph_objs as go
import plotly.figure_factory as ff


old_settings = np.seterr(all='ignore')


class Correlation:
    """
    class to work with matrices
    having sample data
    to compute correlations, plot
    them and make scatter plots
    """

    def __init__(self, matrix_file,
                 corr_method=None,
                 labels=None,
                 remove_outliers=False,
                 skip_zeros=False,
                 log1p=False):

        self.load_matrix(matrix_file)
        self.skip_zeros = skip_zeros
        self.corr_method = corr_method
        self.corr_matrix = None  # correlation matrix
        self.column_order = None
        self.rowCenter = False
        if labels is not None:
            # test that the length of labels
            # corresponds to the length of
            # samples

            self.labels = labels
        self.labels = [toString(x) for x in self.labels]

        if self.matrix.shape[1] == 1:
            # There's nothing that can be done with a single sample
            sys.exit("\nPlease use a matrix with more than one sample\n")

        if skip_zeros is True:
            # remove rows containing only nans or zeros
            # that could be unmappable regions.
            self.remove_rows_of_zeros()

        if remove_outliers is True:
            # remove outliers, otherwise outliers will produce a very
            # high pearson correlation. Unnecessary for spearman correlation
            self.remove_outliers()

        if log1p is True:
            self.matrix = np.log1p(self.matrix)

        if corr_method:
            self.compute_correlation()

    def load_matrix(self, matrix_file):
        """
        loads a matrix file saved using the numpy
        savez method. Two keys are expected:
        'matrix' and 'labels'. The matrix should
        contain one sample per row
        """

        _ma = np.load(matrix_file)
        # matrix:  cols correspond to  samples
        self.matrix = np.asarray(_ma['matrix'].tolist())
        if np.any(np.isnan(self.matrix)):
            num_nam = len(np.flatnonzero(np.isnan(self.matrix.flatten())))
            sys.stderr.write("*Warning*. {} NaN values were found. They will be removed along with the "
                             "corresponding bins in other samples for the computation "
                             "and plotting\n".format(num_nam))

            self.matrix = np.ma.compress_rows(np.ma.masked_invalid(self.matrix))

        self.labels = list(map(toString, _ma['labels']))

        assert len(self.labels) == self.matrix.shape[1], "ERROR, length of labels is not equal " \
                                                         "to length of matrix samples"

    @staticmethod
    def get_outlier_indices(data, max_deviation=200):
        """
        The method is based on the median absolute deviation. See
        Boris Iglewicz and David Hoaglin (1993),
        "Volume 16: How to Detect and Handle Outliers",
        The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.

        returns the list, without the outliers

        The max_deviation=200 is like selecting a z-score
        larger than 200, just that it is based on the median
        and the median absolute deviation instead of the
        mean and the standard deviation.
        """
        median = np.median(data)
        b_value = 1.4826  # value set for a normal distribution
        mad = b_value * np.median(np.abs(data))
        outliers = []
        if mad > 0:
            deviation = abs(data - median) / mad
            """
            outliers = data[deviation > max_deviation]
            print "outliers removed {}".format(len(outliers))
            print outliers
            """
            outliers = np.flatnonzero(deviation > max_deviation)
        return outliers

    def remove_outliers(self, verbose=True):
        """
        get the outliers *per column* using the median absolute
        deviation method

        Returns the filtered matrix
        """

        unfiltered = len(self.matrix)
        to_remove = None
        for col in self.matrix.T:
            outliers = self.get_outlier_indices(col)
            if to_remove is None:
                to_remove = set(outliers)
            else:
                # only set to remove those bins in which
                # the outliers are present in all cases (colums)
                # that's why the intersection is used
                to_remove = to_remove.intersection(outliers)
        if len(to_remove):
            to_keep = [x for x in range(self.matrix.shape[0])
                       if x not in to_remove]
            self.matrix = self.matrix[to_keep, :]
            if verbose:
                sys.stderr.write(
                    "total/filtered/left: "
                    "{}/{}/{}\n".format(unfiltered,
                                        unfiltered - len(to_keep),
                                        len(to_keep)))

        return self.matrix

    def remove_rows_of_zeros(self):
        # remove rows containing all zeros or all nans
        _mat = np.nan_to_num(self.matrix)
        to_keep = _mat.sum(1) != 0

        self.matrix = self.matrix[to_keep, :]

    def save_corr_matrix(self, file_handle):
        """
        saves the correlation matrix
        """
        if self.column_order:
            self.corr_matrix = self.corr_matrix[:, self.column_order][self.column_order]
            self.labels = [self.labels[i] for i in self.column_order]

        self.labels = [toString(x) for x in self.labels]
        file_handle.write("\t'" + "'\t'".join(self.labels) + "'\n")
        fmt = "\t".join(np.repeat('%.4f', self.corr_matrix.shape[1])) + "\n"
        i = 0
        for row in self.corr_matrix:
            file_handle.write(
                "'%s'\t" % self.labels[i] + fmt % tuple(row))
            i += 1

    def compute_correlation(self):
        """
        computes spearman or pearson
        correlation for the samples in the matrix

        The matrix should contain the values of each sample per column
        that's why the transpose is used.

        >>> matrix = np.array([[1, 2, 3, np.nan],
        ...                    [1, 2, 3, 4],
        ...                    [6, 4, 3, 1]]).T
        >>> np.savez_compressed("/tmp/test_matrix.npz", matrix=matrix, labels=['a', 'b', 'c'])

        >>> c = Correlation("/tmp/test_matrix.npz", corr_method='pearson')

        the results should be  as in R

        >>> c.compute_correlation().filled(np.nan)
        array([[ 1.        ,  1.        , -0.98198051],
               [ 1.        ,  1.        , -0.98198051],
               [-0.98198051, -0.98198051,  1.        ]])
        >>> c.corr_method = 'spearman'
        >>> c.corr_matrix = None
        >>> c.compute_correlation()
        array([[ 1.,  1., -1.],
               [ 1.,  1., -1.],
               [-1., -1.,  1.]])
        """
        if self.corr_matrix is not None:
            return self.corr_matrix

        num_samples = len(self.labels)
        # initialize correlation matrix

        if self.corr_method == 'pearson':
            self.corr_matrix = np.ma.corrcoef(self.matrix.T, allow_masked=True)

        else:
            corr_matrix = np.zeros((num_samples, num_samples), dtype='float')
            # do an all vs all correlation using the
            # indices of the upper triangle
            rows, cols = np.triu_indices(num_samples)

            for index in range(len(rows)):
                row = rows[index]
                col = cols[index]
                corr_matrix[row, col] = scipy.stats.spearmanr(self.matrix[:, row], self.matrix[:, col])[0]
            # make the matrix symmetric
            self.corr_matrix = corr_matrix + np.triu(corr_matrix, 1).T

        return self.corr_matrix

    def plotly_correlation(self, corr_matrix, plot_filename, labels, plot_title='',
                           vmax=None, vmin=None, plot_numbers=True,
                           colormap='jet'):
        """plot_correlation, but using plotly"""
        textElement = []
        for row in range(corr_matrix.shape[0]):
            trow = []
            for col in range(corr_matrix.shape[0]):
                if plot_numbers:
                    trow.append("{:0.2f}".format(corr_matrix[row, col]))
                else:
                    trow.append('')
            textElement.append(trow)

        zauto = True
        if vmax is not None or vmin is not None:
            zauto = False

        convertedCmap = convertCmap(colormap)
        fig = ff.create_annotated_heatmap(corr_matrix, x=labels, y=labels, colorscale=convertedCmap, showscale=True, zauto=zauto, zmin=vmin, zmax=vmax, annotation_text=textElement)
        fig.layout['title'] = plot_title
        offline.plot(fig, filename=plot_filename, auto_open=False)

    def plot_correlation(self, plot_filename, plot_title='', vmax=None,
                         vmin=None, colormap='jet', image_format=None,
                         plot_numbers=False, plotWidth=11, plotHeight=9.5):
        """
        plots a correlation using a symmetric heatmap
        """
        num_rows = len(self.labels)
        corr_matrix = self.compute_correlation()
        # set a font size according to figure length
        if num_rows < 6:
            font_size = 14
        elif num_rows > 40:
            font_size = 5
        else:
            font_size = int(14 - 0.25 * num_rows)
        mpl.rcParams.update({'font.size': font_size})
        # set the minimum and maximum values
        if vmax is None:
            vmax = 1
        if vmin is None:
            vmin = 0 if corr_matrix .min() >= 0 else -1

        # Compute and plot dendrogram.
        fig = plt.figure(figsize=(plotWidth, plotHeight))
        plt.suptitle(plot_title)

        axdendro = fig.add_axes([0.02, 0.12, 0.1, 0.66])
        axdendro.set_axis_off()
        y_var = sch.linkage(corr_matrix, method='complete')
        z_var = sch.dendrogram(y_var, orientation='left',
                               link_color_func=lambda k: 'darkred')
        axdendro.set_xticks([])
        axdendro.set_yticks([])
        cmap = plt.get_cmap(colormap)

        # this line simply makes a new cmap, based on the original
        # colormap that goes from 0.0 to 0.9
        # This is done to avoid colors that
        # are too dark at the end of the range that do not offer
        # a good contrast between the correlation numbers that are
        # plotted on black.
        if plot_numbers:
            cmap = pltcolors.LinearSegmentedColormap.from_list(colormap + "clipped",
                                                               cmap(np.linspace(0, 0.9, 10)))

        cmap.set_under((0., 0., 1.))
        # Plot distance matrix.
        axmatrix = fig.add_axes([0.13, 0.1, 0.6, 0.7])
        index = z_var['leaves']
        corr_matrix = corr_matrix[index, :]
        corr_matrix = corr_matrix[:, index]
        if corr_matrix.shape[0] > 30:
            # when there are too many rows it is better to remove
            # the black lines surrounding the boxes in the heatmap
            edge_color = 'none'
        else:
            edge_color = 'black'

        if image_format == "plotly":
            self.plotly_correlation(corr_matrix,
                                    plot_filename,
                                    self.labels,
                                    plot_title=plot_title,
                                    vmax=vmax,
                                    vmin=vmin,
                                    colormap=colormap,
                                    plot_numbers=plot_numbers)
            return

        img_mat = axmatrix.pcolormesh(corr_matrix,
                                      edgecolors=edge_color,
                                      cmap=cmap,
                                      vmax=vmax,
                                      vmin=vmin)
        axmatrix.set_xlim(0, num_rows)
        axmatrix.set_ylim(0, num_rows)

        axmatrix.yaxis.tick_right()
        axmatrix.set_yticks(np.arange(corr_matrix .shape[0]) + 0.5)
        axmatrix.set_yticklabels(np.array(self.labels).astype('str')[index])

        axmatrix.xaxis.set_tick_params(labeltop=True)
        axmatrix.xaxis.set_tick_params(labelbottom=False)
        axmatrix.set_xticks(np.arange(corr_matrix .shape[0]) + 0.5)
        axmatrix.set_xticklabels(np.array(self.labels).astype('str')[index], rotation=45, ha='left')

        axmatrix.tick_params(
            axis='x',
            which='both',
            bottom=False,
            top=False)

        axmatrix.tick_params(
            axis='y',
            which='both',
            left=False,
            right=False)

        # Plot colorbar
        axcolor = fig.add_axes([0.13, 0.065, 0.6, 0.02])
        cobar = plt.colorbar(img_mat, cax=axcolor, orientation='horizontal')
        cobar.solids.set_edgecolor("face")
        if plot_numbers:
            for row in range(num_rows):
                for col in range(num_rows):
                    axmatrix.text(row + 0.5, col + 0.5,
                                  "{:.2f}".format(corr_matrix[row, col]),
                                  ha='center', va='center')

        self.column_order = index
        fig.savefig(plot_filename, format=image_format)
        plt.close()

    def plotly_scatter(self, plot_filename, corr_matrix, plot_title='', minXVal=None, maxXVal=None, minYVal=None, maxYVal=None):
        """Make the scatter plot of a matrix with plotly"""
        n = self.matrix.shape[1]
        self.matrix = self.matrix
        fig = go.Figure()
        domainWidth = 1. / n

        annos = []
        for i in range(n):
            x = domainWidth * (i + 1)
            y = 1 - (domainWidth * i + 0.5 * domainWidth)
            anno = dict(text=self.labels[i], showarrow=False, xref='paper', yref='paper', x=x, y=y, xanchor='right', yanchor='middle')
            annos.append(anno)

        data = []
        zMin = np.inf
        zMax = -np.inf
        for x in range(n):
            xanchor = 'x{}'.format(x + 1)
            base = x * domainWidth
            domain = [base, base + domainWidth]
            if x > 0:
                base = 1 - base
                fig['layout']['xaxis{}'.format(x + 1)] = dict(domain=domain, range=[minXVal, maxXVal], anchor='free', position=base)
            for y in range(0, n):
                yanchor = 'y{}'.format(y + 1)
                if x == 1:
                    base = 1 - y * domainWidth
                    domain = [base - domainWidth, base]
                    fig['layout']['yaxis{}'.format(y + 1)] = dict(domain=domain, range=[minYVal, maxYVal], side='right', anchor='free', position=1.0)

                if x > y:
                    vector1 = self.matrix[:, x]
                    vector2 = self.matrix[:, y]
                    Z, xEdges, yEdges = np.histogram2d(vector1, vector2, bins=50)
                    Z = np.log10(Z)
                    if np.min(Z) < zMin:
                        zMin = np.min(Z)
                    if np.max(Z) > zMax:
                        zMax = np.max(Z)
                    name = '{}={:.2f}'.format(self.corr_method, corr_matrix[x, y])
                    trace = go.Heatmap(z=Z, x=xEdges, y=yEdges, showlegend=False, xaxis=xanchor, yaxis=yanchor, name=name, showscale=False)
                    data.append(trace)

        # Fix the colorbar bounds
        for trace in data:
            trace.update(zmin=zMin, zmax=zMax)
        data[-1]['colorbar'].update(title="log10(instances per bin)", titleside="right")
        data[-1].update(showscale=True)

        fig['data'] = data
        fig['layout'].update(title=plot_title, showlegend=False, annotations=annos)

        offline.plot(fig, filename=plot_filename, auto_open=False)

    def plot_scatter(self, plot_filename, plot_title='', image_format=None, log1p=False, xRange=None, yRange=None):
        """
        Plot the scatter plots of a matrix
        in which each row is a sample
        """

        num_samples = self.matrix.shape[1]
        corr_matrix = self.compute_correlation()
        grids = gridspec.GridSpec(num_samples, num_samples)
        grids.update(wspace=0, hspace=0)
        fig = plt.figure(figsize=(2 * num_samples, 2 * num_samples))
        plt.rcParams['font.size'] = 8.0
        plt.suptitle(plot_title)
        if log1p is True:
            self.matrix = np.log1p(self.matrix)
        min_xvalue = self.matrix.min()
        max_xvalue = self.matrix.max()
        min_yvalue = min_xvalue
        max_yvalue = max_xvalue
        if xRange is not None:
            min_xvalue = xRange[0]
            max_xvalue = xRange[1]
        if yRange is not None:
            min_yvalue = yRange[0]
            max_yvalue = yRange[1]
        if (min_xvalue % 2 == 0 and max_xvalue % 2 == 0) or \
                (min_xvalue % 1 == 0 and max_xvalue % 2 == 1):
            # make one value odd and the other even
            max_xvalue += 1
        if (min_yvalue % 2 == 0 and max_yvalue % 2 == 0) or \
                (min_yvalue % 1 == 0 and max_yvalue % 2 == 1):
            # make one value odd and the other even
            max_yvalue += 1

        # plotly output
        if image_format == 'plotly':
            self.plotly_scatter(plot_filename, corr_matrix, plot_title=plot_title, minXVal=min_xvalue, maxXVal=max_xvalue, minYVal=min_yvalue, maxYVal=max_yvalue)
            return

        rows, cols = np.triu_indices(num_samples)

        for index in range(len(rows)):
            row = rows[index]
            col = cols[index]
            if row == col:
                # add titles as
                # empty plot in the diagonal
                ax = fig.add_subplot(grids[row, col])
                ax.text(0.5, 0.5, self.labels[row],
                        verticalalignment='center',
                        horizontalalignment='center',
                        fontsize=10, fontweight='bold',
                        transform=ax.transAxes)
                ax.set_axis_off()
                continue

            ax = fig.add_subplot(grids[row, col])

            vector1 = self.matrix[:, row]
            vector2 = self.matrix[:, col]

            ax.text(0.2, 0.8, "{}={:.2f}".format(self.corr_method,
                                                 corr_matrix[row, col]),
                    horizontalalignment='left',
                    transform=ax.transAxes)
            ax.get_yaxis().set_tick_params(
                which='both',
                left=False,
                right=False,
                direction='out')

            ax.get_xaxis().set_tick_params(
                which='both',
                top=False,
                bottom=False,
                direction='out')
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_rotation('45')

            if col != num_samples - 1:
                ax.set_yticklabels([])
            else:
                ax.yaxis.tick_right()
                ax.get_yaxis().set_tick_params(
                    which='both',
                    left=False,
                    right=True,
                    direction='out')
            if col - row == 1:
                ax.xaxis.tick_bottom()
                ax.get_xaxis().set_tick_params(
                    which='both',
                    top=False,
                    bottom=True,
                    direction='out')
                for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_rotation('45')

            else:
                ax.set_xticklabels([])

            ax.hist2d(vector2, vector1, bins=200, cmin=0.1)
            ax.set_xlim(min_xvalue, max_xvalue)
            ax.set_ylim(min_yvalue, max_yvalue)

        plt.savefig(plot_filename, format=image_format)
        plt.close()

    def plotly_pca(self, plotFile, Wt, pvar, PCs, eigenvalues, cols, plotTitle):
        """
        A plotly version of plot_pca, that's called by it to do the actual plotting
        """
        fig = go.Figure()
        fig['layout']['xaxis1'] = {'domain': [0.0, 0.48], 'anchor': 'x1', 'title': 'PC{} ({:4.1f}% of var. explained)'.format(PCs[0], 100.0 * pvar[PCs[0] - 1])}
        fig['layout']['yaxis1'] = {'domain': [0.0, 1.0], 'anchor': 'x1', 'title': 'PC{} ({:4.1f}% of var. explained)'.format(PCs[1], 100.0 * pvar[PCs[1] - 1])}
        fig['layout']['xaxis2'] = {'domain': [0.52, 1.0], 'title': 'Principal Component'}
        fig['layout']['yaxis2'] = {'domain': [0.0, 1.0], 'anchor': 'x2', 'title': 'Eigenvalue', 'rangemode': 'tozero', 'showgrid': False}
        fig['layout']['yaxis3'] = {'domain': [0.0, 1.0], 'anchor': 'x2', 'title': 'Cumulative variability', 'rangemode': 'tozero', 'side': 'right', 'overlaying': 'y2'}
        fig['layout'].update(title=plotTitle)

        # PCA
        if cols is not None:
            colors = itertools.cycle(cols)
        n = len(self.labels)
        data = []
        for i in range(n):
            trace = go.Scatter(x=[Wt[PCs[0] - 1, i]],
                               y=[Wt[PCs[1] - 1, i]],
                               mode='marker',
                               xaxis='x1',
                               yaxis='y1',
                               name=self.labels[i])
            trace['marker'].update(size=20)
            if cols is not None:
                trace['marker'].update(color=next(colors))
            data.append(trace)

        # Scree plot
        trace = go.Bar(showlegend=False,
                       name='Eigenvalues',
                       x=range(1, n + 1),
                       y=eigenvalues[:n],
                       xaxis='x2',
                       yaxis='y2')
        data.append(trace)

        # Cumulative variability
        trace = go.Scatter(showlegend=False,
                           x=range(1, n + 1),
                           y=pvar.cumsum()[:n],
                           mode='lines+markers',
                           name='Cumulative variability',
                           xaxis='x2',
                           yaxis='y3',
                           line={'color': 'red'},
                           marker={'symbol': 'circle-open-dot', 'color': 'black'})
        data.append(trace)

        annos = []
        annos.append({'yanchor': 'bottom', 'xref': 'paper', 'xanchor': 'center', 'yref': 'paper', 'text': 'PCA', 'y': 1.0, 'x': 0.25, 'font': {'size': 16}, 'showarrow': False})
        annos.append({'yanchor': 'bottom', 'xref': 'paper', 'xanchor': 'center', 'yref': 'paper', 'text': 'Scree plot', 'y': 1.0, 'x': 0.75, 'font': {'size': 16}, 'showarrow': False})

        fig['data'] = data
        fig['layout']['annotations'] = annos
        offline.plot(fig, filename=plotFile, auto_open=False)

    def plot_pca(self, plot_filename=None, PCs=[1, 2], plot_title='', image_format=None, log1p=False, plotWidth=5, plotHeight=10, cols=None, marks=None):
        """
        Plot the PCA of a matrix

        Returns the matrix of plotted values.
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(plotWidth, plotHeight))

        # Filter
        m = self.matrix
        rvs = m.var(axis=1)
        if self.transpose:
            m = m[np.nonzero(rvs)[0], :]
            rvs = rvs[np.nonzero(rvs)[0]]
        if self.ntop > 0 and m.shape[0] > self.ntop:
            m = m[np.argpartition(rvs, -self.ntop)[-self.ntop:], :]
            rvs = rvs[np.argpartition(rvs, -self.ntop)[-self.ntop:]]

        # log2 (if requested)
        if self.log2:
            self.matrix = np.log2(self.matrix + 0.01)

        # Row center / transpose
        if self.rowCenter and not self.transpose:
            _ = self.matrix.mean(axis=1)
            self.matrix -= _[:, None]
        if self.transpose:
            m = m.T

        # Center and scale
        m2 = (m - np.mean(m, axis=0))
        m2 /= np.std(m2, axis=0, ddof=1)  # Use the unbiased std. dev.

        # SVD
        U, s, Vh = np.linalg.svd(m2, full_matrices=False, compute_uv=True)  # Is full_matrices ever needed?

        # % variance, eigenvalues
        eigenvalues = s**2
        variance = eigenvalues / float(np.max([1, m2.shape[1] - 1]))
        pvar = variance / variance.sum()

        # Weights/projections
        Wt = Vh
        if self.transpose:
            # Use the projected coordinates for the transposed matrix
            Wt = np.dot(m2, Vh.T).T

        if plot_filename is not None:
            n = n_bars = len(self.labels)
            if eigenvalues.size < n:
                n_bars = eigenvalues.size
            markers = itertools.cycle(matplotlib.markers.MarkerStyle.filled_markers)
            if cols is not None:
                colors = itertools.cycle(cols)
            else:
                colors = itertools.cycle(plt.cm.gist_rainbow(np.linspace(0, 1, n)))

            if marks is not None:
                markers = itertools.cycle(marks)

            if image_format == 'plotly':
                self.plotly_pca(plot_filename, Wt, pvar, PCs, eigenvalues, cols, plot_title)
            else:
                ax1.axhline(y=0, color="black", linestyle="dotted", zorder=1)
                ax1.axvline(x=0, color="black", linestyle="dotted", zorder=2)
                for i in range(n):
                    color = next(colors)
                    marker = next(markers)
                    if isinstance(color, np.ndarray):
                        color = pltcolors.to_hex(color, keep_alpha=True)
                    ax1.scatter(Wt[PCs[0] - 1, i], Wt[PCs[1] - 1, i],
                                marker=marker, color=color, s=150, label=self.labels[i], zorder=i + 3)
                if plot_title == '':
                    ax1.set_title('PCA')
                else:
                    ax1.set_title(plot_title)
                ax1.set_xlabel('PC{} ({:4.1f}% of var. explained)'.format(PCs[0], 100.0 * pvar[PCs[0] - 1]))
                ax1.set_ylabel('PC{} ({:4.1f}% of var. explained)'.format(PCs[1], 100.0 * pvar[PCs[1] - 1]))
                lgd = ax1.legend(scatterpoints=1, loc='center left', borderaxespad=0.5,
                                 bbox_to_anchor=(1, 0.5),
                                 prop={'size': 12}, markerscale=0.9)

                # Scree plot
                ind = np.arange(n_bars)  # the x locations for the groups
                width = 0.35        # the width of the bars

                if mpl.__version__ >= "2.0.0":
                    ax2.bar(2 * width + ind, eigenvalues[:n_bars], width * 2)
                else:
                    ax2.bar(width + ind, eigenvalues[:n_bars], width * 2)
                ax2.set_ylabel('Eigenvalue')
                ax2.set_xlabel('Principal Component')
                ax2.set_title('Scree plot')
                ax2.set_xticks(ind + width * 2)
                ax2.set_xticklabels(ind + 1)

                ax3 = ax2.twinx()
                ax3.axhline(y=1, color="black", linestyle="dotted")
                ax3.plot(width * 2 + ind, pvar.cumsum()[:n], "r-")
                ax3.plot(width * 2 + ind, pvar.cumsum()[:n], "wo", markeredgecolor="black")
                ax3.set_ylim([0, 1.05])
                ax3.set_ylabel('Cumulative variability')

                plt.subplots_adjust(top=3.85)
                plt.tight_layout()
                plt.savefig(plot_filename, format=image_format, bbox_extra_artists=(lgd,), bbox_inches='tight')
                plt.close()

        return Wt, eigenvalues
