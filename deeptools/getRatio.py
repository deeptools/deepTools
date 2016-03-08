import numpy as np

old_settings = np.seterr(all='ignore')


def compute_ratio(value1, value2, args):
    value1 = value1 + args['pseudocount']
    value2 = value2 + args['pseudocount']

    ratio = float(value1) / value2
    if args['valueType'] == 'log2':
        ratio = np.log2(ratio)

    elif args['valueType'] == 'reciprocal_ratio':
        # the reciprocal ratio of a/b
        # is a/b if a/b > 1 else -1* b/a
        ratio = ratio if ratio >= 1 else -1.0 / ratio

    return ratio


def getRatio(tileCoverage, args):
    r"""
    The mapreduce method calls this function
    for each tile. The parameters (args) are fixed
    in the main method.

    >>> funcArgs= {'valueType': 'ratio', 'scaleFactors': (1,1), 'pseudocount': 1}
    >>> getRatio([9, 19], funcArgs)
    0.5
    >>> getRatio([0, 0], funcArgs)
    1.0
    >>> getRatio([np.nan, np.nan], funcArgs)
    nan
    >>> getRatio([np.nan, 1.0], funcArgs)
    nan
    >>> funcArgs['valueType'] ='subtract'
    >>> getRatio([20, 10], funcArgs)
    10
    >>> funcArgs['scaleFactors'] = (1, 0.5)
    >>> getRatio([10, 20], funcArgs)
    0.0

    The reciprocal ratio is of a and b is:
    is a/b if a/b > 1 else -1* b/a
    >>> funcArgs['valueType'] ='reciprocal_ratio'
    >>> funcArgs['scaleFactors'] = (1, 1)
    >>> funcArgs['pseudocount'] = 0
    >>> getRatio([2, 1], funcArgs)
    2.0
    >>> getRatio([1, 2], funcArgs)
    -2.0
    >>> getRatio([1, 1], funcArgs)
    1.0
    """

    value1 = args['scaleFactors'][0] * tileCoverage[0]
    value2 = args['scaleFactors'][1] * tileCoverage[1]

    # if any of the two values to compare
    # is nan, return nan
    if np.isnan(value1) or np.isnan(value2):
        return np.nan

    # ratio case
    if args['valueType'] in ['ratio', 'log2', 'reciprocal_ratio']:
        bin_value = compute_ratio(value1, value2, args)

    # non ratio case (diff, sum etc)
    else:
        if args['valueType'] == 'subtract':
            bin_value = value1 - value2
        elif args['valueType'] == 'add':
            bin_value = value1 + value2
        elif args['valueType'] == 'first':
            bin_value = value1
        elif args['valueType'] == 'second':
            bin_value = value2

    return bin_value
