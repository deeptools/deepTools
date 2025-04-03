import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.colors as pltcolors
import textwrap

old_settings = np.seterr(all='ignore')


def plot_single(ax, ma, average_type, color, label, plot_type='lines'):
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
        fill to plot the area between the x axis and the value or any other string to
        just plot the average line.

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
    summary = np.ma.__getattribute__(average_type)(ma, axis=0)
    # only plot the average profiles without error regions
    x = np.arange(len(summary))
    if isinstance(color, np.ndarray):
        color = pltcolors.to_hex(color, keep_alpha=True)
    ax.plot(x, summary, color=color, label=label, alpha=0.9, lw=0.8)
    if plot_type == 'fill':
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

def getProfileTicks(hm, referencePointLabel, startLabel, endLabel, idx):
    """
    returns the position and labelling of the xticks that
    correspond to the heatmap

    As of deepTools 3, the various parameters can be lists, in which case we then need to index things (the idx parameter)

    As of matplotlib 3 the ticks in the heatmap need to have 0.5 added to them.

    As of matplotlib 3.1 there is no longer padding added to all ticks. Reference point ticks will be adjusted by width/2
    or width for spacing and the last half of scaled ticks will be shifed by 1 bin so the ticks are at the beginning of bins.
    """
    w = hm.parameters['bin size']
    b = hm.parameters['upstream']
    a = hm.parameters['downstream']
    if idx is not None:
        w = w[idx]
        b = b[idx]
        a = a[idx]

    try:
        c = hm.parameters['unscaled 5 prime']
        if idx is not None:
            c = c[idx]
    except:
        c = 0
    try:
        d = hm.parameters['unscaled 3 prime']
        if idx is not None:
            d = d[idx]
    except:
        d = 0
    m = hm.parameters['body']
    if idx is not None:
        m = m[idx]

    if b < 1e5:
        quotient = 1000
        symbol = 'Kb'
    else:
        quotient = 1e6
        symbol = 'Mb'

    if m == 0:
        xticks = [(k / w) for k in [0, b - 0.5 * w, b + a - w]]
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

        # We need to subtract the bin size from the last 2 point so they're placed at the beginning of the bin
        xticks_values.append(b + c + m + d - w)
        xtickslabel.append(endLabel)

        if a > 0:
            xticks_values.append(b + c + m + d + a - w)
            xtickslabel.append('{0:.1f}{1}'.format(float(a) / quotient, symbol))

        xticks = [(k / w) for k in xticks_values]
        xticks = [max(x, 0) for x in xticks]

    return xticks, xtickslabel

def justify_text(text, line_width):
    """
    A function to wrap the text with Justify alignment format.
    Ex. If xticks label is longer, then user specify the line width for the text.

    """
    words = text.split()
    lines = textwrap.wrap(text, width=line_width)
    justified_lines = []
    
    for line in lines:
        words_in_line = line.split()
        if len(words_in_line) > 1:
            spaces_needed = line_width - sum(len(w) for w in words_in_line)
            space_between = spaces_needed // (len(words_in_line) - 1)
            extra_spaces = spaces_needed % (len(words_in_line) - 1)
            
            justified_line = ""
            for i, word in enumerate(words_in_line):
                justified_line += word
                if i < len(words_in_line) - 1:
                    justified_line += " " * (space_between + (1 if i < extra_spaces else 0))
            justified_lines.append(justified_line)
        else:
            justified_lines.append(line)
    
    return "\n".join(justified_lines)
