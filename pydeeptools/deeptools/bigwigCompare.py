#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
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
        'per bin, the log2 of the ratio, the sum or the difference.',
        usage='bigwigCompare -b1 sample1.bw -b2 sample2.bw -o log2.bw\n'
        'help: bigwigCompare -h / bigwigCompare --help')

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

    parser.add_argument('--fixedStep',
                        help='Write out all bins (of size --binSize) '
                        'instead of merging neighbouring bins with equal values.',
                        action='store_true')
    return parser


def getType(fname):
    """
    Tries to determine if a file is a wiggle, a bedgraph or a bigWig.
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
        extendPairedEnds=False,
        fixedStep=args.fixedStep)
