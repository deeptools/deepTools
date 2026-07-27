import sys
import itertools
import copy
import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.stats
from deeptools import matplotlib_defaults
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker
import matplotlib.mlab
import matplotlib.markers
import matplotlib.colors as pltcolors
from deeptools.utilities import toString, convertCmap
from scipy.linalg import svd

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
        # Since deeptools 4.0 the labels are encoded and would need to be decoded.
        if _ma['labels'].dtype == np.uint8:
            self.labels = [bytes(x).decode('utf-8').rstrip('\x00') for x in _ma['labels']]
        else:
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

    def plot_correlation(self, plot_filename, plot_title='', vmax=None,
                         vmin=None, colormap='jet', image_format=None,
                         plot_numbers=False, plotWidth=11, plotHeight=9.5, ggplot=False):
        """
        plots a correlation using a symmetric heatmap
        """
        if ggplot:
            plt.style.use('ggplot')
            
        num_rows = len(self.labels)
        corr_matrix = self.compute_correlation()
        # set a font size according to figure length
        if num_rows < 6:
            font_size = 14
        elif num_rows > 40:
            font_size = 5
        else:
            font_size = int(14 - 0.25 * num_rows)
        matplotlib.rcParams.update({'font.size': font_size})
        # set the minimum and maximum values
        if vmax is None:
            vmax = 1
        if vmin is None:
            vmin = 0 if corr_matrix .min() >= 0 else -1

        # Compute and plot dendrogram.
        fig = plt.figure(figsize=(plotWidth, plotHeight))
        plt.suptitle(plot_title)

        axdendro = fig.add_axes([0.015, 0.1, 0.1, 0.7])
        axdendro.set_axis_off()
        y_var = sch.linkage(corr_matrix, method='centroid')
        z_var = sch.dendrogram(y_var, orientation='left',
                               link_color_func=lambda k: 'darkred')
        axdendro.set_xticks([])
        axdendro.set_yticks([])
        cmap = copy.copy(plt.get_cmap(colormap))

        # this line simply makes a new cmap, based on the original
        # colormap that goes from 0.0 to 0.9
        # This is done to avoid colors that
        # are too dark at the end of the range that do not offer
        # a good contrast between the correlation numbers that are
        # plotted on black.
        if plot_numbers:
            cmap = pltcolors.LinearSegmentedColormap.from_list(colormap + "clipped",
                                                               cmap(np.linspace(0, 0.9, 10)))

        cmap = cmap.with_extremes(under=(0., 0., 1.))
        # Plot distance matrix.
        axmatrix = fig.add_axes([0.12, 0.1, 0.6, 0.7])
        index = z_var['leaves']
        corr_matrix = corr_matrix[index, :]
        corr_matrix = corr_matrix[:, index]
        if corr_matrix.shape[0] > 30:
            # when there are too many rows it is better to remove
            # the black lines surrounding the boxes in the heatmap
            edge_color = 'none'
        else:
            edge_color = 'black'

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
        axcolor = fig.add_axes([0.12, 0.065, 0.6, 0.02])
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


    def plot_scatter(self, plot_filename, plot_title='', image_format=None, log1p=False, xRange=None, yRange=None, ggplot=False):
        """
        Plot the scatter plots of a matrix
        in which each row is a sample
        """

        if ggplot:
            plt.style.use('ggplot')

        num_samples = self.matrix.shape[1]
        corr_matrix = self.compute_correlation()
        grids = gridspec.GridSpec(num_samples, num_samples)
        grids.update(wspace=0, hspace=0)
        fig = plt.figure(figsize=(2 * num_samples, 2 * num_samples))
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
            ax.get_xaxis().set_tick_params(
                which='major',
                labelrotation=45)

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
                ax.get_xaxis().set_tick_params(
                    which='major',
                    labelrotation=45)

            else:
                ax.set_xticklabels([])

            ax.set_xlim(min_xvalue, max_xvalue)
            ax.set_ylim(min_yvalue, max_yvalue)
            ax.hist2d(vector2, vector1, bins=200, cmin=0.1)

        plt.savefig(plot_filename, format=image_format)
        plt.close()

    def plot_pca(self, plot_filename=None, PCs=[1, 2], plot_title='', image_format=None, plotWidth=12, plotHeight=10, cols=None, marks=None, add_labels=False, ggplot=False):
        """
        Plot the PCA of a matrix

        Returns the matrix of plotted values.
        """
        if ggplot:
            plt.style.use('ggplot')

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(plotWidth, plotHeight), layout="constrained")

        # Work on a float copy so the transforms below take effect and
        # self.matrix (which callers may reuse) is left untouched.
        m = self.matrix.astype(float, copy=True)

        # log2 (if requested). Applied before variance filtering so the
        # transform actually influences which rows are kept.
        if self.log2:
            m = np.log2(m + 0.01)

        # Row center (incompatible with --transpose, enforced by the CLI).
        if self.rowCenter and not self.transpose:
            m -= m.mean(axis=1)[:, None]

        # Filter to the most variable rows (variance computed post-transform).
        rvs = m.var(axis=1)
        if self.transpose:
            nz = np.nonzero(rvs)[0]
            m = m[nz, :]
            rvs = rvs[nz]
        if self.ntop > 0 and m.shape[0] > self.ntop:
            idx = np.argpartition(rvs, -self.ntop)[-self.ntop:]
            m = m[idx, :]
            rvs = rvs[idx]

        if self.transpose:
            m = m.T

        # Center and scale each column to zero mean and unit variance
        # (equivalent to sklearn's StandardScaler, using population std/ddof=0).
        col_mean = m.mean(axis=0)
        col_std = m.std(axis=0)
        col_std[col_std == 0] = 1.0
        m2 = (m - col_mean) / col_std

        # PCA via SVD of the (re-)centered matrix, mirroring sklearn's PCA().
        n_samples = m2.shape[0]
        X = m2 - m2.mean(axis=0)
        U, S, Vt = svd(X, full_matrices=False)

        # Deterministic sign convention (sklearn's svd_flip, u_based_decision).
        max_abs_cols = np.argmax(np.abs(U), axis=0)
        signs = np.sign(U[max_abs_cols, range(U.shape[1])])
        U *= signs
        Vt *= signs[:, None]

        # Projected coordinates: U * S == X @ V
        Wt = U * S

        # Eigenvalues and % variance explained.
        eigenvalues = (S ** 2) / (n_samples - 1)
        variance = eigenvalues / eigenvalues.sum()
        pvar = variance / variance.sum()

        if self.transpose:
            # With samples as observations, U * S already gives each sample's
            # projection onto the PCs (rows=samples, cols=components). Orient as
            # (components, samples) to match the indexing used below.
            Wt = Wt.T

        if plot_filename is not None:
            n = n_bars = len(self.labels)
            if eigenvalues.size < n:
                n_bars = eigenvalues.size
            # The requested principal components must exist.
            if max(PCs) > eigenvalues.size:
                sys.exit("Cannot plot PC{}: only {} principal component(s) are "
                         "available. Reduce --PCs or increase --ntop.\n".format(max(PCs), eigenvalues.size))
            # In the untransposed layout each point is a sample indexed along
            # the component axis, so there must be at least as many components
            # as samples (i.e. enough usable rows / a large enough --ntop).
            if not self.transpose and Wt.shape[1] < n:
                sys.exit("Not enough principal components ({}) to plot {} "
                         "samples; increase --ntop to at least the sample "
                         "count.\n".format(Wt.shape[1], n))
            markers = itertools.cycle(matplotlib.markers.MarkerStyle.filled_markers)
            if cols is not None:
                colors = itertools.cycle(cols)
            else:
                colors = itertools.cycle(plt.cm.gist_rainbow(np.linspace(0, 1, n)))

            if marks is not None:
                markers = itertools.cycle(marks)

            for i in range(n):
                color = next(colors)
                marker = next(markers)
                if isinstance(color, np.ndarray):
                    color = pltcolors.to_hex(color, keep_alpha=True)
                ax1.scatter(Wt[PCs[0] - 1, i], Wt[PCs[1] - 1, i],
                            marker=marker, color=color, alpha=0.7, label=self.labels[i], zorder=i)
                if add_labels:
                    ax1.text(Wt[PCs[0] - 1, i] * 1.05, Wt[PCs[1] - 1, i] * 1.05, self.labels[i], color="black", fontsize=6, zorder=i)

            # set limits
            xmin = np.min(Wt[PCs[0] - 1, :])
            xmax = np.max(Wt[PCs[0] - 1, :])
            ymin = np.min(Wt[PCs[1] - 1, :])
            ymax = np.max(Wt[PCs[1] - 1, :])

            # symetric limits
            xabs = max(abs(xmin), abs(xmax))
            yabs = max(abs(ymin), abs(ymax))
            xmin = -xabs
            xmax = xabs
            ymin = -yabs
            ymax = yabs

            # add some space
            xmin -= 0.05 * (xmax - xmin)
            xmax += 0.05 * (xmax - xmin)
            ymin -= 0.05 * (ymax - ymin)
            ymax += 0.05 * (ymax - ymin)
            ax1.set_xlim([xmin, xmax])
            ax1.set_ylim([ymin, ymax])

            # labels
            if plot_title == '':
                ax1.set_title('PCA')
            else:
                ax1.set_title(plot_title)
            ax1.set_xlabel('PC{} ({:4.1f}% of var. explained)'.format(PCs[0], 100.0 * pvar[PCs[0] - 1]))
            ax1.set_ylabel('PC{} ({:4.1f}% of var. explained)'.format(PCs[1], 100.0 * pvar[PCs[1] - 1]))

            if not add_labels:
                if n < 30:
                    ncols = 1
                else:
                    ncols = 2
                lgd = ax1.legend(scatterpoints=1, loc='upper left', borderaxespad=0.,
                                 bbox_to_anchor=(1.05, 1), ncols=ncols,
                                 prop={'size': 8}, markerscale=0.9)
            else:
                lgd = None

            # Scree plot
            ind = np.arange(n_bars) + 1 # the x locations for the groups

            ax2.plot(ind, pvar, "bo-")
            ax2.set_xlabel('Principal Component')
            ax2.set_title('Scree plot')
            ax2.set_xticks(ind)
            ax2.plot(ind, pvar.cumsum()[:n], "ro-")
            ax2.set_ylabel('Variability')

            ax2.axhline(y=1, color="black", linestyle="dotted")

            ax2.legend(['individual', 'accumulative'], loc="center right", labelcolor=["blue", "red"], ncols=1)

            if lgd is not None:
                plt.savefig(plot_filename, format=image_format, bbox_extra_artists=(lgd,))
            else:
                plt.savefig(plot_filename, format=image_format)
            plt.close()

        return Wt, eigenvalues
