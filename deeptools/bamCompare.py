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
        'containing the bin location and the resulting comparison value. By '
        'default, if reads are paired, the fragment length reported in the BAM '
        'file is used. Each mate, however, '
        'is treated independently to avoid a bias when a mixture of concordant '
        'and discordant pairs is present. This means that *each end* will '
        'be extended to match the fragment length.',

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
                          'Default "readCount".',
                          choices=['readCount', 'SES'],
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

    optional.add_argument('--ratio',
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


def get_scale_factors(args):
    if args.ratio == 'subtract':
        # We need raw counts in this case
        normalizeTo1x = args.normalizeTo1x
        normalizeUsingRPKM = args.normalizeUsingRPKM
        args.normalizeTo1x = False
        args.normalizeUsingRPKM = False

    # This is only used if we subtract
    mapped_reads = [None, None]

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
        args.bam = args.bamfile1
        args.scaleFactor = 1.0
        bam1_mapped, _ = get_num_kept_reads(args)
        args.bam = args.bamfile2
        bam2_mapped, _ = get_num_kept_reads(args)
        scale_factors = float(min(bam1_mapped, bam2_mapped)) / np.array([bam1_mapped, bam2_mapped])
        mapped_reads = [bam1_mapped, bam2_mapped]
        if args.verbose:
            print("Size factors using total number "
                  "of mapped reads: {}".format(scale_factors))

    # in case the subtract method is used, the final difference
    # would be normalized according to the given method
    if args.ratio == 'subtract':
        # The next lines identify which of the samples is not scaled down.
        # The normalization using RPKM or normalize to 1x would use
        # as reference such sample. Since the other sample would be
        # scaled to match the un-scaled one, the normalization factor due to RPKM or normalize1x
        # for both samples should be based on the unscaled one.
        # For example, if sample A is unscaled and sample B is scaled by 0.5,
        # then normalizing factor for A to report RPKM read counts
        # is also applied to B.

        if args.scaleFactors is None:
            # check which of the two samples is not scaled down
            if scale_factors[0] == 1:
                args.bam = args.bamfile1
                mapped_reads = mapped_reads[0]
            else:
                args.bam = args.bamfile2
                mapped_reads = mapped_reads[1]
            if mapped_reads is None:
                mapped_reads, _ = get_num_kept_reads(args)

        # Replace the arguments
        args.normalizeTo1x = normalizeTo1x
        args.normalizeUsingRPKM = normalizeUsingRPKM

        if args.scaleFactors is None:
            if args.normalizeTo1x:
                # try to guess fragment length if the bam file contains paired end reads
                from deeptools.getFragmentAndReadSize import get_read_and_fragment_length
                frag_len_dict, read_len_dict = get_read_and_fragment_length(args.bam,
                                                                            return_lengths=False,
                                                                            blackListFileName=args.blackListFileName,
                                                                            numberOfProcessors=args.numberOfProcessors,
                                                                            verbose=args.verbose)
                if args.extendReads:
                    if args.extendReads is True:
                        # try to guess fragment length if the bam file contains paired end reads
                        if frag_len_dict:
                            fragment_length = frag_len_dict['median']
                        else:
                            exit("*ERROR*: library is not paired-end. Please provide an extension length.")
                        if args.verbose:
                            print(("Fragment length based on paired en data "
                                  "estimated to be {}".format(frag_len_dict['median'])))

                    elif args.extendReads < 1:
                        exit("*ERROR*: read extension must be bigger than one. Value give: {} ".format(args.extendReads))
                    elif args.extendReads > 2000:
                        exit("*ERROR*: read extension must be smaller that 2000. Value give: {} ".format(args.extendReads))
                    else:
                        fragment_length = args.extendReads

                else:
                    # set as fragment length the read length
                    fragment_length = int(read_len_dict['median'])
                    if args.verbose:
                        print("Estimated read length is {}".format(int(read_len_dict['median'])))

                current_coverage = float(mapped_reads * fragment_length) / args.normalizeTo1x
                # the coverage scale factor is 1 / coverage,
                coverage_scale_factor = 1.0 / current_coverage
                scale_factors = np.array(scale_factors) * coverage_scale_factor
                if args.verbose:
                    print("Estimated current coverage {}".format(current_coverage))
                    print("Scale factor to convert "
                          "current coverage to 1: {}".format(coverage_scale_factor))
            else:
                # by default normalize using RPKM
                # the RPKM is:
                # Num reads per tile/(total reads (in millions)*tile length in Kb)
                millionReadsMapped = float(mapped_reads) / 1e6
                tileLengthInKb = float(args.binSize) / 1000
                coverage_scale_factor = 1.0 / (millionReadsMapped * tileLengthInKb)
                scale_factors = np.array(scale_factors) * coverage_scale_factor
                if args.verbose:
                    print("Scale factor for RPKM is {0}".format(coverage_scale_factor))

    return scale_factors


def main(args=None):
    """
    The algorithm is composed of two parts.

    1. Using the SES or read counts method, appropriate scaling
       factors are determined to account for sequencing depth differences.

    2. The genome is transversed, scaling the BAM files, and computing
       the log ratio/ratio/difference for bins of fixed width
       given by the user.

    """
    args = process_args(args)

    scale_factors = get_scale_factors(args)
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
                                     verbose=args.verbose
                                     )

    wr.run(FUNC, func_args, args.outFileName, blackListFileName=args.blackListFileName, format=args.outFileFormat, smoothLength=args.smoothLength)

if __name__ == "__main__":
    main()
