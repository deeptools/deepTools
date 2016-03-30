#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse  # to parse command line arguments
from deeptools import parserCommon
from deeptools.getRatio import getRatio
from deeptools import writeBedGraph_bam_and_bw

debug = 0


def parse_arguments(args=None):
    parentParser = parserCommon.getParentArgParse()
    outputParser = parserCommon.output()
    parser = argparse.ArgumentParser(
        parents=[parentParser, outputParser],
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
                        help='small number to avoid x/0. Only useful '
                        'when ratio = log2 or ratio',
                        default=1,
                        type=float,
                        required=False)

    parser.add_argument('--ratio',
                        help='The default is to output the log2ratio of the '
                        'two samples. The reciprocal ratio returns the '
                        'the negative of the inverse of the ratio '
                        'if the ratio is less than 0. The resulting '
                        'values are interpreted as negative fold changes. '
                        '*NOTE*: Only with --ratio subtract can --normalizeTo1x or '
                        '--normalizeUsingRPKM be used. Instead of performing a '
                        'computation using both files, the scaled signal can '
                        'alternatively be output for the first or second file using '
                        'the \'--ratio first\' or \'--ratio second\'',
                        default='log2',
                        choices=['log2', 'ratio', 'subtract', 'add',
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


def main(args=None):
    args = parse_arguments().parse_args(args)

    if args.scaleFactors:
        scaleFactors = [float(x) for x in args.scaleFactors.split(":")]
    else:
        scaleFactors = [1, 1]

    # the getRatio function is called and receives
    # the function_args per each tile that is considered
    FUNC = getRatio
    function_args = {'valueType': args.ratio,
                     'scaleFactors': scaleFactors,
                     'pseudocount': args.pseudocount}

    writeBedGraph_bam_and_bw.writeBedGraph(
        [(args.bigwig1, 'bigwig'),
         (args.bigwig2, 'bigwig')],
        args.outFileName, 0, FUNC,
        function_args, tileSize=args.binSize, region=args.region,
        blackListFileName=args.blackListFileName,
        numberOfProcessors=args.numberOfProcessors,
        format=args.outFileFormat,
        smoothLength=False,
        missingDataAsZero=not args.skipNonCoveredRegions,
        extendPairedEnds=False)
