#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse  # to parse command line arguments
import sys
import multiprocessing
import os
import numpy as np
from deeptools import parserCommon
from deeptools import writeBedGraph_bam_and_bw
import deeptools.deepBlue as db

debug = 0


def parse_arguments(args=None):
    parentParser = parserCommon.getParentArgParse()
    outputParser = parserCommon.output()
    dbParser = parserCommon.deepBlueOptionalArgs()
    parser = argparse.ArgumentParser(
        parents=[parentParser, outputParser, dbParser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='This tool average multiple bigWig files based on the number '
        'of mapped reads. To average the bigWig files, the genome is '
        'partitioned into bins of equal size, then the scores '
        'in each bigwig file are computed per bin.'
        'These scores are averaged and scaleFactors can be applied before the average.',
        usage='bigwigAverage -b sample1.bw sample2.bw -o outfile.bw\n'
        'help: bigwigAverage -h / bigwigAverage --help')

    # define the arguments
    parser.add_argument('--bigwigs', '-b',
                        metavar='Bigwig files',
                        help='Bigwig files separated by space.',
                        nargs='+',
                        required=True)

    parser.add_argument('--scaleFactors',
                        help='Set this parameter to multipy the bigwig values '
                        'by a constant. The format is '
                        'scaleFactor1:scaleFactor2:scaleFactor3 etc. '
                        'For example 0.7:1 to scale the first bigwig file '
                        'by 0.7 while not scaling the second bigwig file',
                        default=None,
                        required=False)

    parser.add_argument('--skipNonCoveredRegions', '--skipNAs',
                        help='This parameter determines if non-covered regions (regions without a score) '
                        'in the bigWig files should be skipped. The default is to treat those '
                        'regions as having a value of zero. '
                        'The decision to skip non-covered regions '
                        'depends on the interpretation of the data. Non-covered regions '
                        'in a bigWig file may represent repetitive regions that should '
                        'be skipped. Alternatively, the interpretation of non-covered regions as '
                        'zeros may be wrong and this option should be used ',
                        action='store_true')

    return parser


def getType(fname):
    """
    Tries to determine if a file is a wiggle file from deepBlue or a bigWig file.
    Returns 'wiggle' if the file name ends with .wig, otherwise 'bigwig'
    """
    if fname.endswith(".wig") or fname.endswith(".wiggle"):
        return "wiggle"
    elif fname.lower().endswith(".bedgraph") or fname.endswith(".bdg"):
        return "bedgraph"
    else:
        return "bigwig"


def average(tileCoverage, args):
    r"""
    The mapreduce method calls this function
    for each tile. The parameters (args) are fixed
    in the main method.

    >>> funcArgs= {'scaleFactors': (1,1)}
    >>> average([1, 2], funcArgs)
    1.5
    >>> funcArgs= {'scaleFactors': (1,0.5)}
    >>> average([1, 2], funcArgs)
    1.0
    >>> funcArgs= {'scaleFactors': (1,0.5,0.1,0.2)}
    >>> average([1, 2, 3, 12], funcArgs)
    1.175
    >>> average([1, 2, 3, np.nan], funcArgs)
    nan
    """

    norm_values = [args['scaleFactors'][i] * cov for i, cov in enumerate(tileCoverage)]

    return np.mean(norm_values)


def main(args=None):
    args = parse_arguments().parse_args(args)
    if len(sys.argv) == 1:
        parse_arguments().print_help()
        sys.exit()

    nFiles = len(args.bigwigs)

    if args.scaleFactors:
        scaleFactors = [float(x) for x in args.scaleFactors.split(":")]
        if len(scaleFactors) == 1:
            scaleFactors = scaleFactors * nFiles
        elif len(scaleFactors) != nFiles:
            raise argparse.ArgumentTypeError(
                "Format of scaleFactors is factor or factor1:factor2... as many as bigwig files. "
                "There are {} bigwigs and {} factors."
                "The value given ( {} ) is not valid".format(nFiles, len(scaleFactors), args.scaleFactors))
    else:
        scaleFactors = [1] * nFiles

    # the average function is called and receives
    # the function_args per each tile that is considered
    FUNC = average
    function_args = {'scaleFactors': scaleFactors}

    # Preload deepBlue files, which need to then be deleted
    deepBlueFiles = []
    for idx, fname in enumerate(args.bigwigs):
        if db.isDeepBlue(fname):
            deepBlueFiles.append([fname, idx])
    if len(deepBlueFiles) > 0:
        sys.stderr.write("Preloading the following deepBlue files: {}\n".format(",".join([x[0] for x in deepBlueFiles])))
        foo = db.deepBlue(deepBlueFiles[0][0], url=args.deepBlueURL, userKey=args.userKey)
        regs = db.makeChromTiles(foo)
        for x in deepBlueFiles:
            x.extend([args, regs])
        if len(deepBlueFiles) > 1 and args.numberOfProcessors > 1:
            pool = multiprocessing.Pool(args.numberOfProcessors)
            res = pool.map_async(db.preloadWrapper, deepBlueFiles).get(9999999)
        else:
            res = list(map(db.preloadWrapper, deepBlueFiles))

        # substitute the file names with the temp files
        for (ftuple, r) in zip(deepBlueFiles, res):
            args.bigwigs[ftuple[1]] = r
        deepBlueFiles = [[x[0], x[1]] for x in deepBlueFiles]
        del regs

    writeBedGraph_bam_and_bw.writeBedGraph(
        [(b, getType(b)) for b in args.bigwigs],
        args.outFileName, 0, FUNC,
        function_args, tileSize=args.binSize, region=args.region,
        blackListFileName=args.blackListFileName,
        verbose=args.verbose,
        numberOfProcessors=args.numberOfProcessors,
        skipZeroOverZero=False,
        format=args.outFileFormat,
        smoothLength=False,
        missingDataAsZero=not args.skipNonCoveredRegions,
        extendPairedEnds=False)

    # Clean up temporary bigWig files, if applicable
    if not args.deepBlueKeepTemp:
        for k, v in deepBlueFiles:
            os.remove(args.bigwigs[v])
    else:
        for k, v in deepBlueFiles:
            foo = args.bigwigs[v]
            print("{} is stored in {}".format(k, foo))
