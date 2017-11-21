#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse  # to parse command line arguments
import numpy as np

# my packages
from deeptools import writeBedGraph
from deeptools.SES_scaleFactor import estimateScaleFactor
from deeptools import parserCommon
from deeptools import bamHandler
from deeptools.getRatio import getRatio
from deeptools.getScaleFactor import get_num_kept_reads
from deeptools.getScaleFactor import get_scale_factor
debug = 0
old_settings = np.seterr(all='ignore')


def parseArguments():
    parentParser = parserCommon.getParentArgParse()
    bamParser = parserCommon.read_options()
    normalizationParser = parserCommon.normalization_options()
    requiredArgs = getRequiredArgs()
    optionalArgs = getOptionalArgs()
    outputParser = parserCommon.output()
    parser = argparse.ArgumentParser(
        parents=[requiredArgs, outputParser, optionalArgs,
                 parentParser, normalizationParser, bamParser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='This tool compares two BAM files based on the number of '
        'mapped reads. To compare the BAM files, the genome is partitioned '
        'into bins of equal size, then the number of reads found in each bin'
        ' is counted per file, and finally a summary value is '
        'reported. This value can be the ratio of the number of reads per '
        'bin, the log2 of the ratio, or the difference. This tool can '
        'normalize the number of reads in each BAM file using the SES method '
        'proposed by Diaz et al. (2012) "Normalization, bias correction, and '
        'peak calling for ChIP-seq". Statistical Applications in Genetics '
        'and Molecular Biology, 11(3). Normalization based on read counts '
        'is also available. The output is either a bedgraph or bigWig file '
        'containing the bin location and the resulting comparison value. '
        'Note that *each end* in a pair (for paired-end reads) is treated '
        'independently. If this is undesirable, then use the --samFlagInclude '
        'or --samFlagExclude options.',

        usage=' bamCompare -b1 treatment.bam -b2 control.bam -o log2ratio.bw',

        add_help=False)

    return parser


def getRequiredArgs():
    parser = argparse.ArgumentParser(add_help=False)

    required = parser.add_argument_group('Required arguments')

    # define the arguments
    required.add_argument('--bamfile1', '-b1',
                          metavar='BAM file',
                          help='Sorted BAM file 1. Usually the BAM file '
                          'for the treatment.',
                          required=True)

    required.add_argument('--bamfile2', '-b2',
                          metavar='BAM file',
                          help='Sorted BAM file 2. Usually the BAM '
                          'file for the control.',
                          required=True)

    return parser


def getOptionalArgs():

    parser = argparse.ArgumentParser(add_help=False)
    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")

    optional.add_argument('--scaleFactorsMethod',
                          help='Method to use to scale the samples. '
                          'If a method is specified, then it will be used to compensate '
                          'for sequencing depth differences between the samples. '
                          'As an alternative, this can be set to None and an option from '
                          '--normalizeUsing <method> can be used.',
                          choices=['readCount', 'SES', 'None'],
                          default='readCount')

    optional.add_argument('--sampleLength', '-l',
                          help='*Only relevant when SES is chosen for the '
                          'scaleFactorsMethod.* To compute the SES, specify '
                          'the length (in bases) of the regions (see --numberOfSamples) '
                          'that will be randomly sampled to calculate the scaling factors. '
                          'If you do not have a good sequencing depth for '
                          'your samples consider increasing the sampling '
                          'regions\' size to minimize the probability '
                          'that zero-coverage regions are used.',
                          default=1000,
                          type=int)

    optional.add_argument('--numberOfSamples', '-n',
                          help='*Only relevant when SES is chosen for the '
                          'scaleFactorsMethod.* Number of samplings taken '
                          'from the genome to compute the scaling factors.',
                          default=1e5,
                          type=int)

    optional.add_argument('--scaleFactors',
                          help='Set this parameter manually to avoid the computation of '
                          'scaleFactors. The format is scaleFactor1:scaleFactor2.'
                          'For example, --scaleFactor 0.7:1 will cause the first BAM file to'
                          'be multiplied by 0.7, while not scaling '
                          'the second BAM file (multiplication with 1).',
                          default=None,
                          required=False)

    optional.add_argument('--ratio', '--method',
                          help='The default is to output the log2 ratio of the '
                          'two samples. The reciprocal ratio returns the '
                          'the negative of the inverse of the ratio '
                          'if the ratio is less than 0. The resulting '
                          'values are interpreted as negative fold changes. '
                          'Instead of performing a computation using both files, the scaled signal can '
                          'alternatively be output for the first or second file using '
                          'the \'--ratio first\' or \'--ratio second\'',
                          default='log2',
                          choices=['log2', 'ratio', 'subtract', 'add', 'mean',
                                   'reciprocal_ratio', 'first', 'second'],
                          required=False)

    optional.add_argument('--pseudocount',
                          help='small number to avoid x/0. Only useful '
                          'together with --ratio log2 or --ratio ratio .',
                          default=1,
                          type=float,
                          required=False)

    return parser


def process_args(args=None):
    args = parseArguments().parse_args(args)

    if args.smoothLength and args.smoothLength <= args.binSize:
        print("Warning: the smooth length given ({}) is smaller than the bin "
              "size ({}).\n\n No smoothing will be "
              "done".format(args.smoothLength,
                            args.binSize))
        args.smoothLength = None

    if not args.ignoreForNormalization:
        args.ignoreForNormalization = []
    return args

# get_scale_factors function is used for scaling in bamCompare
# while get_scale_factor is used for depth normalization


def get_scale_factors(args):

    if args.scaleFactors:
        scale_factors = list(map(float, args.scaleFactors.split(":")))
    elif args.scaleFactorsMethod == 'SES':
        scalefactors_dict = estimateScaleFactor(
            [args.bamfile1, args.bamfile2],
            args.sampleLength, args.numberOfSamples,
            1,
            blackListFileName=args.blackListFileName,
            numberOfProcessors=args.numberOfProcessors,
            verbose=args.verbose,
            chrsToSkip=args.ignoreForNormalization)

        scale_factors = scalefactors_dict['size_factors']

        if args.verbose:
            bam1 = bamHandler.openBam(args.bamfile1)
            bam2 = bamHandler.openBam(args.bamfile2)

            print("Size factors using SES: {}".format(scale_factors))
            print("%s regions of size %s where used " %
                  (scalefactors_dict['sites_sampled'],
                   args.sampleLength))

            print("ignoring filtering/blacklists, size factors if the number of mapped "
                  "reads would have been used:")
            print(tuple(
                float(min(bam1.mapped, bam2.mapped)) / np.array([bam1.mapped, bam2.mapped])))
            bam1.close()
            bam2.close()

    elif args.scaleFactorsMethod == 'readCount':
        # change the scaleFactor to 1.0
        args.scaleFactor = 1.0
        # get num of kept reads for bam file 1
        args.bam = args.bamfile1
        bam1_mapped, _ = get_num_kept_reads(args)
        # get num of kept reads for bam file 2
        args.bam = args.bamfile2
        bam2_mapped, _ = get_num_kept_reads(args)

        mapped_reads = [bam1_mapped, bam2_mapped]

        # new scale_factors (relative to min of two bams)
        scale_factors = float(min(bam1_mapped, bam2_mapped)) / np.array(mapped_reads)
        if args.verbose:
            print("Size factors using total number "
                  "of mapped reads: {}".format(scale_factors))

    elif args.scaleFactorsMethod == 'None':
        scale_factors = None

    return scale_factors


def main(args=None):
    """
    The algorithm is composed of two steps.


    1. Per-sample scaling / depth Normalization:
     + If scaling is used (using the SES or read counts method), appropriate scaling
       factors are determined to account for sequencing depth differences.
     + Optionally scaling can be turned off and individual samples could be depth normalized using
       RPKM, BPM or CPM methods

    2. Ratio calculation between two bam files:
     + The genome is transversed and computing
       the log ratio/ratio/difference etc. for bins of fixed width
       given by the user.

    """
    args = process_args(args)

    if args.normalizeTo1x is not None:
        print("Warning! RPGC normalization (--normalizeTo1x) is not supported with bamCompare. Ignored..")
        args.normalizeTo1x = None

    scale_factors = get_scale_factors(args)
    if scale_factors is None:
        # check whether one of the depth norm methods are selected
        if args.normalizeUsing is not None:
            args.scaleFactor = 1.0
            # if a normalization is required then compute the scale factors
            args.bam = args.bamfile1
            scale_factor_bam1 = get_scale_factor(args)
            args.bam = args.bamfile2
            scale_factor_bam2 = get_scale_factor(args)
            scale_factors = [scale_factor_bam1, scale_factor_bam2]
        else:
            scale_factors = [1, 1]

    if args.verbose:
        print("Individual scale factors are {0}".format(scale_factors))

    # the getRatio function is called and receives
    # the func_args per each tile that is considered
    FUNC = getRatio
    func_args = {'valueType': args.ratio,
                 'scaleFactors': scale_factors,
                 'pseudocount': args.pseudocount
                 }

    wr = writeBedGraph.WriteBedGraph([args.bamfile1, args.bamfile2], args.binSize, 0,
                                     stepSize=args.binSize,
                                     region=args.region,
                                     numberOfProcessors=args.numberOfProcessors,
                                     extendReads=args.extendReads,
                                     blackListFileName=args.blackListFileName,
                                     minMappingQuality=args.minMappingQuality,
                                     ignoreDuplicates=args.ignoreDuplicates,
                                     center_read=args.centerReads,
                                     zerosToNans=args.skipNonCoveredRegions,
                                     samFlag_include=args.samFlagInclude,
                                     samFlag_exclude=args.samFlagExclude,
                                     minFragmentLength=args.minFragmentLength,
                                     maxFragmentLength=args.maxFragmentLength,
                                     chrsToSkip=args.ignoreForNormalization,
                                     verbose=args.verbose
                                     )

    wr.run(FUNC, func_args, args.outFileName, blackListFileName=args.blackListFileName, format=args.outFileFormat, smoothLength=args.smoothLength)


if __name__ == "__main__":
    main()
