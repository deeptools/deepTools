import numpy as np


def getRatio(tileCoverage, args):
    r"""
    The mapreduce method calls this function
    for each tile. The parameters (args) are fixed
    in the main method.

    >>> funcArgs= {'missingDataAsZero': True, 'valueType': 'ratio',
    ... 'scaleFactors': (1,1), 'p1': 0.1, 'p2': 0.1}
    >>> getRatio([10,20], funcArgs)
    0.5
    >>> getRatio([0,0], funcArgs)
    1.0
    >>> getRatio([np.nan,np.nan], funcArgs)
    1.0
    >>> getRatio([np.nan,1.0], funcArgs)
    0.09090909090909091
    >>> funcArgs['missingDataAsZero'] = False
    >>> getRatio([10,np.nan], funcArgs)
    nan
    >>> funcArgs['valueType'] ='subtract'
    >>> getRatio([20,10], funcArgs)
    10
    >>> funcArgs['scaleFactors'] = (1, 0.5)
    >>> getRatio([10,20], funcArgs)
    0.0
    >>> funcArgs['valueType'] ='reciprocal_ratio'
    >>> funcArgs['scaleFactors'] = (1, 1)
    >>> getRatio([20,10], funcArgs)
    2.0
    >>> getRatio([10,20], funcArgs)
    -2.0
    """
    if not args['missingDataAsZero']:
        if np.isnan(args['scaleFactors'][0]) or \
                np.isnan(args['scaleFactors'][1]):
            return np.nan

    value1 = args['scaleFactors'][0] * tileCoverage[0]
    value2 = args['scaleFactors'][1] * tileCoverage[1]

    if args['missingDataAsZero'] is True:
        if np.isnan(value1):
            value1 = 0
        if np.isnan(value2):
            value2 = 0

    # case when both tile coverage counts are zero
    if (value1 == 0.0 or np.isnan(value1)) and \
            (value2 == 0.0 or np.isnan(value2)):
        if args['missingDataAsZero'] is True:
            if args['valueType'] == 'subtract' or \
                    args['valueType'] == 'log2' or \
                    args['valueType'] == 'add':
                ratio = 0.0
            else:
                ratio = 1.00
        else:
            ratio = np.nan
    else:
        if args['valueType'] == 'subtract':
            ratio = value1 - value2
        elif args['valueType'] == 'add':
            ratio = value1 + value2
        else:
            # the pseudocount is only useful when ratios are considered
            if value2 == 0.0 or np.isnan(value2) or \
                    value1 == 0.0 or np.isnan(value1):
                ratio = float(value1 + args['p1']) / (value2 + args['p2'])
            else:
                ratio = float(value1) / value2

            if args['valueType'] == 'log2':
                if ratio == 0:
                    ratio = np.nan
                else:
                    ratio = np.log2(ratio)
            elif args['valueType'] == 'reciprocal_ratio':
                # the reciprocal ratio of a/b
                # is a/b if a/b > 1 else -1* b/a
                ratio = ratio if ratio > 1 else -1.0 / ratio

    return ratio
