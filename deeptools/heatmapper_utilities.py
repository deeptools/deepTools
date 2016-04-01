import numpy as np
import matplotlib.colors as pltcolors

old_settings = np.seterr(all='ignore')


def plot_single(ax, ma, average_type, color, label, plot_type='simple'):
    """
    Adds a line to the plot in the given ax using the specified method

    Parameters
    ----------
    ax : matplotlib axis
        matplotlib axis
    ma : numpy array
        numpy array The data on this matrix is summarized according
        to the `average_type` argument.
    average_type : str
        string values are sum mean median min max std
    color : str
        a valid color: either a html color name, hex
        (e.g #002233), RGB + alpha tuple or list or RGB tuple or list
    label : str
        label
    plot_type: str
        type of plot. Either 'se' for standard error, 'std' for
        standard deviation, 'overlapped_lines' to plot each line of the matrix,
        fill to plot the area between the x axis and the value or None, just to
        plot the average line.

    Returns
    -------
    ax
        matplotlib axis

    Examples
    --------

    >>> import matplotlib.pyplot as plt
    >>> import os
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> matrix = np.array([[1,2,3],
    ...                    [4,5,6],
    ...                    [7,8,9]])
    >>> ax = plot_single(ax, matrix -2, 'mean', color=[0.6, 0.8, 0.9], label='fill light blue', plot_type='fill')
    >>> ax = plot_single(ax, matrix, 'mean', color='blue', label='red')
    >>> ax = plot_single(ax, matrix + 5, 'mean', color='red', label='red', plot_type='std')
    >>> ax = plot_single(ax, matrix + 10, 'mean', color='#cccccc', label='gray se', plot_type='se')
    >>> ax = plot_single(ax, matrix + 20, 'mean', color=(0.9, 0.5, 0.9), label='violet', plot_type='std')
    >>> ax = plot_single(ax, matrix + 30, 'mean', color=(0.9, 0.5, 0.9, 0.5), label='violet with alpha', plot_type='std')
    >>> leg = ax.legend()
    >>> plt.savefig("/tmp/test.pdf")
    >>> plt.close()
    >>> fig = plt.figure()
    >>> os.remove("/tmp/test.pdf")


    """
    summary = np.__getattribute__(average_type)(ma, axis=0)
    # only plot the average profiles without error regions
    x = np.arange(len(summary))
    ax.plot(x, summary, color=color, label=label, alpha=0.9)
    if plot_type == 'fill':
        pass
        ax.fill_between(x, summary, facecolor=color, alpha=0.6, edgecolor='none')

    if plot_type in ['se', 'std']:
        if plot_type == 'se':  # standard error
            std = np.std(ma, axis=0) / np.sqrt(ma.shape[0])
        else:
            std = np.std(ma, axis=0)

        alpha = 0.2
        # an alpha channel has to be added to the color to fill the area
        # between the mean (or median etc.) and the std or se
        f_color = pltcolors.colorConverter.to_rgba(color, alpha)

        ax.fill_between(x, summary, summary + std, facecolor=f_color, edgecolor='none')
        ax.fill_between(x, summary, summary - std, facecolor=f_color, edgecolor='none')

    ax.set_xlim(0, max(x))

    return ax


def getProfileTicks(hm, referencePointLabel, startLabel, endLabel):
    """
    returns the position and labelling of the xticks that
    correspond to the heatmap
    """
    w = hm.parameters['bin size']
    b = hm.parameters['upstream']
    a = hm.parameters['downstream']
    try:
        c = hm.parameters['unscaled 5 prime']
    except:
        c = 0
    try:
        d = hm.parameters['unscaled 3 prime']
    except:
        d = 0
    m = hm.parameters['body']
    tickPlotAdj = 0.5

    if b < 1e5:
        quotient = 1000
        symbol = 'Kb'
    else:
        quotient = 1e6
        symbol = 'Mb'

    if m == 0:
        xticks = [(k / w) - tickPlotAdj for k in [0, b, b + a]]
        xtickslabel = ['{0:.1f}'.format(-(float(b) / quotient)),
                       referencePointLabel,
                       '{0:.1f}{1}'.format(float(a) / quotient, symbol)]

    else:
        xticks_values = [0]
        xtickslabel = []

        # only if upstream region is set, add a x tick
        if b > 0:
            xticks_values.append(b)
            xtickslabel.append('{0:.1f}'.format(-(float(b) / quotient)))

        xtickslabel.append(startLabel)

        # set the x tick for the body parameter, regardless if
        # upstream is 0 (not set)
        if c > 0:
            xticks_values.append(b + c)
            xtickslabel.append("")

        if d > 0:
            xticks_values.append(b + c + m)
            xtickslabel.append("")

        xticks_values.append(b + c + m + d)
        xtickslabel.append(endLabel)

        if a > 0:
            xticks_values.append(b + c + m + d + a)
            xtickslabel.append('{0:.1f}{1}'.format(float(a) / quotient, symbol))

        xticks = [(k / w) - tickPlotAdj for k in xticks_values]

    return xticks, xtickslabel
