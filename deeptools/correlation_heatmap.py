from matplotlib import use as mplt_use
mplt_use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.cluster.hierarchy as sch
from matplotlib import rcParams

old_settings = np.seterr(all='ignore')


def plot_correlation(corr_matrix, labels, plotFileName, vmax=None,
                     vmin=None, colormap='jet', image_format=None,
                     plot_numbers=False, plot_title=''):

    num_rows = corr_matrix.shape[0]

    # set a font size according to figure length
    if num_rows < 6:
        font_size = 14
    elif num_rows > 40:
        font_size = 5
    else:
        font_size = int(14 - 0.25 * num_rows)
    rcParams.update({'font.size': font_size})
    # set the minimum and maximum values
    if vmax is None:
        vmax = 1
    if vmin is None:
        vmin = 0 if corr_matrix.min() >= 0 else -1

    # Compute and plot dendrogram.
    fig = plt.figure(figsize=(11, 9.5))
    if plot_title:
        plt.suptitle(plot_title)
    axdendro = fig.add_axes([0.02, 0.12, 0.1, 0.66])
    axdendro.set_axis_off()
    y_var = sch.linkage(corr_matrix, method='complete')
    z_var = sch.dendrogram(y_var, orientation='right',
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
    img_mat = axmatrix.pcolormesh(corr_matrix,
                                  edgecolors='black',
                                  cmap=cmap,
                                  vmax=vmax,
                                  vmin=vmin)
    axmatrix.set_xlim(0, num_rows)
    axmatrix.set_ylim(0, num_rows)

    axmatrix.yaxis.tick_right()
    axmatrix.set_yticks(np.arange(corr_matrix.shape[0]) + 0.5)
    axmatrix.set_yticklabels(np.array(labels).astype('str')[index])

#    axmatrix.xaxis.set_label_position('top')
    axmatrix.xaxis.set_tick_params(labeltop='on')
    axmatrix.xaxis.set_tick_params(labelbottom='off')
    axmatrix.set_xticks(np.arange(corr_matrix.shape[0]) + 0.5)
    axmatrix.set_xticklabels(np.array(labels).astype('str')[index],
                             rotation=45,
                             ha='left')

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

    #    axmatrix.set_xticks([])
    # Plot colorbar.
    axcolor = fig.add_axes([0.13, 0.065, 0.6, 0.02])
    cobar = plt.colorbar(img_mat, cax=axcolor, orientation='horizontal')
    cobar.solids.set_edgecolor("face")
    if plot_numbers:
        for row in range(num_rows):
            for col in range(num_rows):
                axmatrix.text(row + 0.5, col + 0.5,
                              "{:.2f}".format(corr_matrix[row, col]),
                              ha='center', va='center')

    fig.savefig(plotFileName, format=image_format)
    fig.close()
