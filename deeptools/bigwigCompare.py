#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse  # to parse command line arguments
import sys
import multiprocessing
import os
from deeptools import parserCommon
from deeptools.getRatio import getRatio
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
        description='This tool compares two bigWig files based on the number '
        'of mapped reads. To compare the bigWig files, the genome is '
        'partitioned into bins of equal size, then the number of reads found '
        'in each BAM file are counted per bin and finally a summary '
        'value is reported. This value can be the ratio of the number of reads'
        'per bin, the log2 of the ratio, the sum or the difference.')

    # define the arguments
    parser.add_argument('--bigwig1', '-b1',
                        metavar='Bigwig file',
                        help='Bigwig file 1. Usually the file for the '
                        'treatment.',
                        required=True)

    parser.add_argument('--bigwig2', '-b2',
                        metavar='Bigwig file',
                        help='Bigwig file 2. Usually the file for the '
                        'control.',
                        required=True)

    parser.add_argument('--scaleFactors',
                        help='Set this parameter to multipy the bigwig values '
                        'by a constant. The format is '
                        'scaleFactor1:scaleFactor2. '
                        'For example 0.7:1 to scale the first bigwig file '
                        'by 0.7 while not scaling the second bigwig file',
                        default=None,
                        required=False)

    parser.add_argument('--pseudocount',
                        help='A small number to avoid x/0. Only useful '
                        'together with --operation log2 or --operation ratio. '
                        'You can specify different values as pseudocounts for '
                        'the numerator and the denominator by providing two '
                        'values (the first value is used as the numerator '
                        'pseudocount and the second the denominator pseudocount). (Default: %(default)s)',
                        default=1,
                        nargs='+',
                        action=parserCommon.requiredLength(1, 2),
                        type=float,
                        required=False)

    parser.add_argument('--skipZeroOverZero',
                        help='Skip bins where BOTH BAM files lack coverage. '
                        'This is determined BEFORE any applicable pseudocount '
                        'is added.',
                        action='store_true')

    parser.add_argument('--operation',
                        help='The default is to output the log2ratio of the '
                        'two samples. The reciprocal ratio returns the '
                        'the negative of the inverse of the ratio '
                        'if the ratio is less than 0. The resulting '
                        'values are interpreted as negative fold changes. '
                        'Instead of performing a '
                        'computation using both files, the scaled signal can '
                        'alternatively be output for the first or second file using '
                        'the \'--operation first\' or \'--operation second\' (Default: %(default)s)',
                        default='log2',
                        choices=['log2', 'ratio', 'subtract', 'add', 'mean',
                                 'reciprocal_ratio', 'first', 'second'],
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
    elif fname.endswith(".bedgraph"):
        return "bedgraph"
    else:
        return "bigwig"


def main(args=None):
    args = parse_arguments().parse_args(args)

    if args.scaleFactors:
        scaleFactors = [float(x) for x in args.scaleFactors.split(":")]
    else:
        scaleFactors = [1, 1]

    if not isinstance(args.pseudocount, list):
        args.pseudocount = [args.pseudocount]

    if len(args.pseudocount) == 1:
        args.pseudocount *= 2

    # the getRatio function is called and receives
    # the function_args per each tile that is considered
    FUNC = getRatio
    function_args = {'valueType': args.operation,
                     'scaleFactors': scaleFactors,
                     'pseudocount': args.pseudocount}

    # Preload deepBlue files, which need to then be deleted
    deepBlueFiles = []
    for idx, fname in enumerate([args.bigwig1, args.bigwig2]):
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
            if ftuple[1] == 0:
                args.bigwig1 = r
            else:
                args.bigwig2 = r
        deepBlueFiles = [[x[0], x[1]] for x in deepBlueFiles]
        del regs

    writeBedGraph_bam_and_bw.writeBedGraph(
        [(args.bigwig1, getType(args.bigwig1)),
         (args.bigwig2, getType(args.bigwig2))],
        args.outFileName, 0, FUNC,
        function_args, tileSize=args.binSize, region=args.region,
        blackListFileName=args.blackListFileName,
        verbose=args.verbose,
        numberOfProcessors=args.numberOfProcessors,
        skipZeroOverZero=args.skipZeroOverZero,
        format=args.outFileFormat,
        smoothLength=False,
        missingDataAsZero=not args.skipNonCoveredRegions,
        extendPairedEnds=False)

    # Clean up temporary bigWig files, if applicable
    if not args.deepBlueKeepTemp:
        for k, v in deepBlueFiles:
            if v == 0:
                os.remove(args.bigwig1)
            else:
                os.remove(args.bigwig2)
    else:
        for k, v in deepBlueFiles:
            foo = args.bigwig1
            if v == 1:
                foo = args.bigwig2
            print("{} is stored in {}".format(k, foo))
