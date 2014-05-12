import numpy as np

def plot_single(ax, ma, average_type, color, label,
                 plot_type='simple'):

    """
    Adds a line to the plot using the specified metod
    """

    sumry = np.__getattribute__(average_type)(ma, axis=0)
    # only plot the average profiles without error regions
    if plot_type != 'overlapped_lines':
        ax.plot(sumry, color=color, label=label, alpha=0.9)
        x = np.arange(len(sumry))
    if plot_type == 'fill':
        ax.fill_between(x, sumry, facecolor=color, alpha=0.6)

    elif plot_type == 'se':  #standard error
        std = np.std(ma, axis=0)/np.sqrt(ma.shape[0])
        alpha = 0.2
        if type(color) == type((0,0)):  # check of color is tuple
            # add the alpha channed to the color tuple
            f_color = [c for c in color[:-1]] + [alpha]
        else:
            f_color = pltcolors.colorConverter.to_rgba(color, alpha)
        # ideally the edgecolor should be None,
        # but when generating a pdf image an edge is still
        # drawn.
        ax.fill_between(x, sumry, sumry + std, facecolor=f_color,
                        edgecolor=f_color, lw=0.01)
        ax.fill_between(x, sumry, sumry - std, facecolor=f_color,
                        edgecolor=f_color, lw=0.01)

    elif plot_type == 'overlapped_lines':
        ax.patch.set_facecolor('black')
        for row in ma:
            ax.plot(row, 'yellow', alpha=0.1)
        x = np.arange(len(row))
    ax.set_xlim(0, max(x))


def getProfileTicks(hm, referencePointLabel, startLabel, endLabel):
    """
    returns the position and labelling of the xticks that
    correspond to the heatmap
    """
    w = hm.parameters['bin size']
    b = hm.parameters['upstream']
    a = hm.parameters['downstream']
    m = hm.parameters['body']
    tickPlotAdj = 0.5

    if b < 1e5:
        quotient = 1000
        symbol = 'Kb'
    if b >= 1e5:
        quotient = 1e6
        symbol = 'Mb'

    if m == 0:
        xticks = [(k / w) - tickPlotAdj for k in [0, b, b  + a]]
        xtickslabel = ['{0:.1f}'.format(-(float(b) / quotient)),
                       referencePointLabel,
                       '{0:.1f}{1}'.format(float(a) / quotient, symbol)]

    else:
        xticks_values = [0]
        xtickslabel = []

        # only if upstream region is set, add a x tick
        if hm.parameters['upstream'] > 0:
            xticks_values.append(b)
            xtickslabel.append('{0:.1f}'.format(-(float(b) / quotient)))

        # set the x tick for the body parameter, regardless if
        # upstream is 0 (not set)
        xticks_values.append(b + m)
        xtickslabel.append(startLabel)
        xtickslabel.append(endLabel)
        if a > 0:
            xticks_values.append(b + m + a)
            xtickslabel.append('{0:.1f}{1}'.format(float(a) / quotient, symbol))

        xticks = [(k / w) - tickPlotAdj for k in xticks_values]

    return xticks, xtickslabel
