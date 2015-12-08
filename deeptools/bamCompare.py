#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys
import argparse  # to parse command line arguments
import numpy as np

# my packages
from deeptools import writeBedGraph
from deeptools.SES_scaleFactor import estimateScaleFactor
from deeptools import parserCommon
from deeptools import bamHandler
from deeptools.getRatio import getRatio

debug = 0


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
        'into bins of equal size, then the number of reads found in each BAM '
        'file is counted for such bins and finally a summarizing value is '
        'reported. This value can be the ratio of the number of reads per '
        'bin, the log2 of the ratio or the difference. This tool can '
        'normalize the number of reads on each BAM file using the SES method '
        'proposed by Diaz et al. (2012). "Normalization, bias correction, and '
        'peak calling for ChIP-seq". Statistical applications in genetics '
        'and molecular biology, 11(3). Normalization based on read counts '
        'is also available. The output is either a bedgraph or a bigwig file '
        'containing the bin location and the resulting comparison values. By '
        'default, if reads are mated, the fragment length reported in the BAM '
        'file is used. In the case of paired-end mapping each read mate '
        'is treated independently to avoid a bias when a mixture of concordant '
        'and discordant pairs is present. This means that *each end* will '
        'be extended to match the fragment length.',

        usage='An example usage is:\n %(prog)s '
        '-b1 treatment.bam -b2 control.bam -o log2ratio.bw',

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

    optional.add_argument('--bamIndex1', '-bai1',
                          help='Index for the bam file1 . Default is '
                          'to consider a the path of the bam file adding '
                          'the .bai suffix.',
                          metavar='bam file index')

    optional.add_argument('--bamIndex2', '-bai2',
                          help='Index for the bam file1. Default is to '
                          'consider the path of the bam file adding the .bai '
                          'suffix.',
                          metavar='bam file index')

    optional.add_argument('--scaleFactorsMethod',
                          help='Method to use to scale the samples. '
                          'Default "readCount".',
                          choices=['readCount', 'SES'],
                          default='readCount')

    optional.add_argument('--sampleLength', '-l',
                          help='*Only relevant when SES is chosen for the '
                          'scaleFactorsMethod.* To compute the SES, specify '
                          'the length (in bp) of the regions (see --numberOfSamples) '
                          'which will be randomly sampled to calculate the scaling factors. '
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
                          help='Set this parameter to avoid the computation of '
                          'scaleFactors. The format is '
                          'scaleFactor1:scaleFactor2. For example 0.7:1 to '
                          'scale the first BAM file by 0.7 while not scaling '
                          'the second BAM file',
                          default=None,
                          required=False)

    optional.add_argument('--ratio',
                          help='The default is to output the log2ratio between the '
                          'two samples. The reciprocal ratio returns the '
                          'the negative of the inverse of the ratio '
                          'if the ratio is less than 0. The resulting '
                          'values are interpreted as negative fold changes. '
                          '*NOTE*: Only when --ratio subtract, the options --normalizeTo1x or '
                          '--normalizeUsingRPKM can be used.',
                          default='log2',
                          choices=['log2', 'ratio', 'subtract', 'add',
                                   'reciprocal_ratio'],
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
    args.missingDataAsZero = True if args.missingDataAsZero == 'yes' else False

    if args.smoothLength and args.smoothLength <= args.binSize:
        print "Warning: the smooth length given ({}) is smaller than the bin "\
            "size ({}).\n\n No smoothing will be "\
            "done".format(args.smoothLength,
                          args.binSize)
        args.smoothLength = None

    if not args.ignoreForNormalization:
        args.ignoreForNormalization = []
    return args


def get_scale_factors(args):

    bam1 = bamHandler.openBam(args.bamfile1, args.bamIndex1)
    bam2 = bamHandler.openBam(args.bamfile2, args.bamIndex2)

    bam1_mapped = parserCommon.bam_total_reads(bam1, args.ignoreForNormalization)
    bam2_mapped = parserCommon.bam_total_reads(bam2, args.ignoreForNormalization)

    if args.scaleFactors:
        scale_factors = map(float, args.scaleFactors.split(":"))
    else:
        if args.scaleFactorsMethod == 'SES':
            scalefactors_dict = estimateScaleFactor(
                [bam1.filename, bam2.filename],
                args.sampleLength, args.numberOfSamples,
                1,
                numberOfProcessors=args.numberOfProcessors,
                verbose=args.verbose,
                chrsToSkip=args.ignoreForNormalization)

            scale_factors = scalefactors_dict['size_factors']

            if args.verbose:
                print "Size factors using SES: {}".format(scale_factors)
                print "%s regions of size %s where used " % \
                    (scalefactors_dict['sites_sampled'],
                     args.sampleLength)

                print "size factor if the number of mapped " \
                    "reads would have been used:"
                print tuple(
                    float(min(bam1.mapped, bam2.mapped)) / np.array([bam1.mapped, bam2.mapped]))

        elif args.scaleFactorsMethod == 'readCount':
            scale_factors = \
                float(min(bam1_mapped, bam2_mapped)) / np.array([bam1_mapped, bam2_mapped])
            if args.verbose:
                print "Size factors using total number " \
                    "of mapped reads: {}".format(scale_factors)

    # in case the subtract method is used, the final difference
    # would be normalized according to the given method
    if args.ratio == 'subtract':
        # The next lines identify which of the samples is not scaled down.
        # The normalization using RPKM or normalize to 1x would use
        # as reference such sample. Since the other sample would be
        # scaled to match the un-scaled one, the normalization factor
        # for both samples should be based on the unscaled one.
        # For example, if sample A is unscaled and sample B is scaled by 0.5,
        # then normalizing factor for A to report RPKM read counts
        # is also applied to B.
        if scale_factors[0] == 1:
            mappedReads = bam1_mapped
            bamfile = args.bamfile1
            bamindex = args.bamIndex1
        else:
            mappedReads = bam2_mapped
            bamfile = args.bamfile2
            bamindex = args.bamIndex2

        if args.scaleFactors is None:
            if args.normalizeTo1x:
                from deeptools.getFragmentAndReadSize import get_read_and_fragment_length
                frag_len_dict, read_len_dict = get_read_and_fragment_length(bamfile, bamindex,
                                                                            return_lengths=False,
                                                                            numberOfProcessors=args.numberOfProcessors,
                                                                            verbose=args.verbose)
                if args.fragmentLength:
                    if frag_len_dict['mean'] != 0 and \
                                    abs(args.fragmentLength - frag_len_dict['median']) > frag_len_dict['std']:
                        sys.stderr.write("*Warning*:\nThe fragment length provided ({}) does not match the fragment "
                                         "length estimated from the bam file: {}\n".format(args.fragmentLength,
                                                                                         int(frag_len_dict['median'])))

                    fragment_length = args.fragmentLength

                else:
                    # set as fragment length the read length
                    if args.verbose:
                        print "Estimated read length is {}".format(int(read_len_dict['median']))
                    fragment_length = int(read_len_dict['median'])


                current_coverage = \
                                 float(mappedReads * fragment_length) / args.normalizeTo1x
                # the coverage scale factor is 1 / coverage,
                coverage_scale_factor = 1.0 / current_coverage
                scale_factors = np.array(scale_factors) * coverage_scale_factor
                if args.verbose:
                    print "Estimated current coverage {}".format(current_coverage)
                    print "Scale factor to convert " \
                          "current coverage to 1: {}".format(coverage_scale_factor)
            else:
                # by default normalize using RPKM
                # the RPKM is:
                # Num reads per tile/(total reads (in millions)*tile length in Kb)
                millionReadsMapped = float(mappedReads)  / 1e6
                tileLengthInKb = float(args.binSize) / 1000
                coverage_scale_factor = 1.0 / (millionReadsMapped * tileLengthInKb)
                scale_factors = np.array(scale_factors) * coverage_scale_factor

                if args.verbose:
                    print "scale factor for   "
                    "RPKM is {0}".format(coverage_scale_factor)

    return scale_factors


def main(args=None):
    """
    The algorithm is composed of two parts.

    1. Using the SES or mapped reads method, appropriate scaling
       factors are determined.

    2. The genome is transversed, scaling the BAM files, and computing
       the log ratio/ratio/difference for bins of fixed width
       given by the user.

    """
    args = process_args(args)

    scale_factors = get_scale_factors(args)
    if args.verbose:
        print "Individual scale factors are {0}".format(scale_factors)

    # the getRatio function is called and receives
    # the func_args per each tile that is considered
    FUNC = getRatio
    func_args = {'missingDataAsZero': args.missingDataAsZero,
                 'valueType': args.ratio,
                 'scaleFactors': scale_factors,
                 'pseudocount': args.pseudocount
                 }

    wr = writeBedGraph.WriteBedGraph([args.bamfile1, args.bamfile2], args.binSize, 0,
                                     stepSize=args.binSize,
                                     region=args.region,
                                     numberOfProcessors=args.numberOfProcessors,
                                     extendReads=args.extendReads,
                                     minMappingQuality=args.minMappingQuality,
                                     ignoreDuplicates=args.ignoreDuplicates,
                                     center_read=args.centerReads,
                                     zerosToNans=True,
                                     samFlag_include=args.samFlagInclude,
                                     samFlag_exclude=args.samFlagExclude,
                                     verbose=args.verbose
                                     )

    wr.run(FUNC, func_args,  args.outFileName, format=args.outFileFormat, smooth_length=args.smoothLength)

if __name__ == "__main__":
    main()
