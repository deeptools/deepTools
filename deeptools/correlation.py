import sys
import itertools
import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker
import matplotlib.mlab
import matplotlib.markers

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
        if labels is not None:
            # test that the length of labels
            # corresponds to the length of
            # samples

            self.labels = labels

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

        self.labels = _ma['labels']

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
            self.corr_matrix[self.column_order, self.column_order]
            self.labels = self.labels[self.column_order]

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

    def plot_correlation(self, plot_fiilename, plot_title='', vmax=None,
                         vmin=None, colormap='jet', image_format=None,
                         plot_numbers=False):
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
        fig = plt.figure(figsize=(11, 9.5))
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
            cmap = cmap.from_list(colormap + "clipped",
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

        axmatrix.xaxis.set_tick_params(labeltop='on')
        axmatrix.xaxis.set_tick_params(labelbottom='off')
        axmatrix.set_xticks(np.arange(corr_matrix .shape[0]) + 0.5)
        axmatrix.set_xticklabels(np.array(self.labels).astype('str')[index], rotation=45, ha='left')

        axmatrix.tick_params(
            axis='x',
            which='both',
            bottom='off',
            top='off')

        axmatrix.tick_params(
            axis='y',
            which='both',
            left='off',
            right='off')

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
        fig.savefig(plot_fiilename, format=image_format)
        plt.close()

    def plot_scatter(self, plot_fiilename, plot_title='', image_format=None, log1p=False):
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
        min_value = self.matrix.min()
        max_value = self.matrix.max()
        if (min_value % 2 == 0 and max_value % 2 == 0) or \
                (min_value % 1 == 0 and max_value % 2 == 1):
            # make one value odd and the other even
            max_value += 1

        if log1p:
            major_locator = matplotlib.ticker.FixedLocator(list(range(min_value, max_value, 2)))
            minor_locator = matplotlib.ticker.FixedLocator(list(range(min_value, max_value, 1)))

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
            if log1p:
                ax.xaxis.set_major_locator(major_locator)
                ax.xaxis.set_minor_locator(minor_locator)
                ax.yaxis.set_major_locator(major_locator)
                ax.yaxis.set_minor_locator(minor_locator)

            vector1 = self.matrix[:, row]
            vector2 = self.matrix[:, col]

            ax.text(0.2, 0.8, "{}={:.2f}".format(self.corr_method,
                                                 corr_matrix[row, col]),
                    horizontalalignment='left',
                    transform=ax.transAxes)
            ax.get_yaxis().set_tick_params(
                which='both',
                left='off',
                right='off',
                direction='out')

            ax.get_xaxis().set_tick_params(
                which='both',
                top='off',
                bottom='off',
                direction='out')
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_rotation('45')

            if col != num_samples - 1:
                ax.set_yticklabels([])
            else:
                ax.yaxis.tick_right()
                ax.get_yaxis().set_tick_params(
                    which='both',
                    left='off',
                    right='on',
                    direction='out')
            if col - row == 1:
                ax.xaxis.tick_bottom()
                ax.get_xaxis().set_tick_params(
                    which='both',
                    top='off',
                    bottom='on',
                    direction='out')
                for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_rotation('45')

            else:
                ax.set_xticklabels([])

            ax.hist2d(vector1, vector2, bins=200, cmin=0.1)
            # downsample for plotting
    #        choice_idx = np.random.randint(0, len(vector1),min(len(vector1), 500000))
    #        ax.plot(vector1[choice_idx], vector2[choice_idx], '.', markersize=1,
    #                    alpha=0.3, color='darkblue',
    #                    markeredgecolor=None)

    #        ax.set_ylim(min_value, max_value)
    #        ax.set_xlim(min_value,max_value)
            ax.set_ylim(min_value, ax.get_ylim()[1])
            ax.set_xlim(min_value, ax.get_xlim()[1])

        plt.savefig(plot_fiilename, format=image_format)
        plt.close()

    def plot_pca(self, plot_filename, plot_title='', image_format=None, log1p=False):
        """
        Plot the PCA of a matrix
        """

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 10))
        # PCA
        mlab_pca = matplotlib.mlab.PCA(self.matrix)
        n = len(self.labels)
        markers = itertools.cycle(matplotlib.markers.MarkerStyle.filled_markers)
        colors = itertools.cycle(plt.cm.gist_rainbow(np.linspace(0, 1, n)))

        ax1.axhline(y=0, color="black", linestyle="dotted", zorder=1)
        ax1.axvline(x=0, color="black", linestyle="dotted", zorder=2)
        for i in range(n):
            ax1.scatter(mlab_pca.Wt[0, i], mlab_pca.Wt[1, i],
                        marker=next(markers), color=next(colors), s=150, label=self.labels[i], zorder=i + 3)
        if plot_title == '':
            ax1.set_title('PCA')
        else:
            ax1.set_title(plot_title)
        ax1.set_xlabel('PC1')
        ax1.set_ylabel('PC2')
        lgd = ax1.legend(scatterpoints=1, loc='center left', borderaxespad=0.5,
                         bbox_to_anchor=(1, 0.5),
                         prop={'size': 12}, markerscale=0.9)

        # Scree plot
        eigenvalues = mlab_pca.s

        cumulative = []
        c = 0
        for x in mlab_pca.fracs:
            c += x
            cumulative.append(c)

        ind = np.arange(n)  # the x locations for the groups
        width = 0.35        # the width of the bars

        ax2.bar(width + ind, eigenvalues, width * 2)
        ax2.set_ylabel('Eigenvalue')
        ax2.set_xlabel('Factors')
        ax2.set_title('Scree plot')
        ax2.set_xticks(ind + width * 2)
        ax2.set_xticklabels(ind + 1)

        ax3 = ax2.twinx()
        ax3.axhline(y=1, color="black", linestyle="dotted")
        ax3.plot(width * 2 + ind, cumulative[0:], "r-")
        ax3.plot(width * 2 + ind, cumulative[0:], "wo")
        ax3.set_ylim([0, 1.05])
        ax3.set_ylabel('Cumulative variability')

        plt.subplots_adjust(top=3.85)
        plt.tight_layout()
        plt.savefig(plot_filename, format=image_format, bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close()
