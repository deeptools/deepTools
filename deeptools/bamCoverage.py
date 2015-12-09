#!/usr/bin/env python
# -*- coding: utf-8 -*-

# own tools
import argparse
import sys
from deeptools import writeBedGraph  # This should be made directly into a bigWig
from deeptools import parserCommon
from deeptools import bamHandler

debug = 0


def parseArguments():
    parentParser = parserCommon.getParentArgParse()
    bamParser = parserCommon.read_options()
    requiredArgs = get_required_args()
    optionalArgs = get_optional_args()
    outputParser = parserCommon.output()
    parser = \
        argparse.ArgumentParser(
            parents=[requiredArgs, outputParser, optionalArgs,
                     parentParser, bamParser],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            description='Given a BAM file, this tool generates a bigWig or '
            'bedGraph file of read or fragment coverage. The method first '
            'calculates the number of reads (either extended to match the '
            'fragment length or not) that overlap each bin in the genome.\n'
            'The resulting read counts can be '
            'normalized using either a given scaling factor or the RPKM formula, '
            'or to get a 1x depth of coverage (RPGC).\n'
            'In the case of paired-end mapping, each read mate is treated '
            'independently to avoid a bias when a mixture of concordant and '
            'discordant pairs is present. This means that *each mate* will '
            'be extended to match the fragment length.',
            usage='An example usage is: %(prog)s -b signal.bam -o signal.bw',
            add_help=False)

    return parser


def get_required_args():
    parser = argparse.ArgumentParser(add_help=False)

    required = parser.add_argument_group('Required arguments')

    # define the arguments
    required.add_argument('--bam', '-b',
                          help='BAM file to process',
                          metavar='BAM file',
                          required=True)

    return parser


def get_optional_args():

    parser = argparse.ArgumentParser(add_help=False)
    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")

    optional.add_argument('--bamIndex', '-bai',
                          help='Index for the BAM file. Default is to consider '
                          'the path of the BAM file adding the .bai suffix.',
                          metavar='BAM file index')

    optional.add_argument('--scaleFactor',
                          help='Indicate a number that you would like to use. When used in combination '
                          'with --normalizeTo1x or --normalizeUsingRPKM, the computed '
                          'scaling factor will be multiplied by the given scale factor.',
                          default=1.0,
                          type=float,
                          required=False)

    optional.add_argument('--normalizeTo1x',
                          help='Report read coverage normalized to 1x '
                          'sequencing depth (also known as Reads Per Genomic '
                          'Content (RPGC)). Sequencing depth is defined as: '
                          '(total number of mapped reads * fragment length) / '
                          'effective genome size.\nThe scaling factor used '
                          'is the inverse of the sequencing depth computed '
                          'for the sample to match the 1x coverage. '
                          'To use this option, the '
                          'effective genome size has to be indicated after the '
                          'option. The effective genome size is the portion '
                          'of the genome that is mappable. Large fractions of '
                          'the genome are stretches of NNNN that should be '
                          'discarded. Also, if repetitive regions were not '
                          'included in the mapping of reads, the effective '
                          'genome size needs to be adjusted accordingly. '
                          'Common values are: mm9: 2,150,570,000; '
                          'hg19:2,451,960,000; dm3:121,400,000 and ce10:93,260,000. '
                          'See Table 2 of http://www.plosone.org/article/info:doi/10.1371/journal.pone.0030377 '
                          'or http://www.nature.com/nbt/journal/v27/n1/fig_tab/nbt.1518_T1.html '
                          'for several effective genome sizes.',
                          metavar='EFFECTIVE GENOME SIZE LENGTH',
                          default=None,
                          type=int,
                          required=False)

    optional.add_argument('--normalizeUsingRPKM',
                          help='Use Reads Per Kilobase per Million reads to '
                          'normalize the number of reads per bin. The formula '
                          'is: RPKM (per bin) =  number of reads per bin / '
                          '( number of mapped reads (in millions) * bin '
                          'length (kb) ). Each read is considered independently,'
                          'if you want to only count either of the mate pairs in'
                          'paired-end data, use the --samFlag option.',
                          action='store_true',
                          required=False)

    optional.add_argument('--ignoreForNormalization', '-ignore',
                          help='A list of chromosome names separated by spaces '
                          'containing those chromosomes that should be excluded '
                          'for computing the normalization. This is useful when considering '
                          'samples with unequal coverage across chromosomes like male '
                          'samples. An usage examples is  --ignoreForNormalization chrX chrM.',
                          nargs='+')

    optional.add_argument('--missingDataAsZero',
                          default="yes",
                          choices=["yes", "no"],
                          help='If set to "yes", missing data will be treated as zero. '
                          'If set to "no", missing data will be ignored '
                          'and not included in the output file. Missing '
                          'data is defined as those bins for which '
                          'no overlapping reads are found.')

    optional.add_argument('--smoothLength',
                          metavar="INT bp",
                          help='The smooth length defines a window, larger than '
                          'the binSize, to average the number of reads. For '
                          'example, if the --binSize is set to 20 bp and the '
                          '--smoothLength is set to 60 bp, then, for each '
                          'bin, the average of the bin and its left and right '
                          'neighbors is considered. Any value smaller than the '
                          '--binSize will be ignored and no smoothing will be '
                          'applied.',
                          type=int)

    optional.add_argument('--MNase',
                          help='Determine nucleosome positions from MNase-seq data. '
                          'Only 3 nucleotides at the center of each fragment are counted. '
                          'The fragment ends are defined by the two mate reads. '
                          '*NOTE*: Requires paired-end data.',
                          action='store_true')

    return parser


def scaleFactor(string):
    try:
        scaleFactor1, scaleFactor2 = string.split(":")
        scaleFactors = (float(scaleFactor1), float(scaleFactor2))
    except:
        raise argparse.ArgumentTypeError(
            "Format of scaleFactors is factor1:factor2. "
            "The value given ( {} ) is not valid".format(string))
    return scaleFactors


def process_args(args=None):
    args = parseArguments().parse_args(args)

    if args.scaleFactor != 1:
        args.normalizeTo1x = None
    if args.smoothLength and args.smoothLength <= args.binSize:
        print "Warning: the smooth length given ({}) is smaller than the bin "\
            "size ({}).\n\n No smoothing will "\
            "be done".format(args.smoothLength,
                             args.binSize)
        args.smoothLength = None

    return args


def get_scale_factor(args):

    scale_factor = args.scaleFactor
    bamHandle = bamHandler.openBam(args.bam, args.bamIndex)
    bam_mapped = parserCommon.bam_total_reads(bamHandle, args.ignoreForNormalization)

    if args.normalizeTo1x:
        # try to guess fragment length if the bam file contains paired end reads
        from deeptools.getFragmentAndReadSize import get_read_and_fragment_length
        frag_len_dict, read_len_dict = get_read_and_fragment_length(args.bam, args.bamIndex,
                                                                    return_lengths=False,
                                                                    numberOfProcessors=args.numberOfProcessors,
                                                                    verbose=args.verbose)
        if args.extendReads and args.extendReads is True:
            fragment_length = args.frag_len_dict['median']

        elif args.extendReads:
            if frag_len_dict['mean'] != 0 and abs(args.extendReads - frag_len_dict['median']) > frag_len_dict['std']:
                sys.stderr.write("*Warning*:\nThe extend reads length provided ({}) does not match the fragment "
                                 "length estimated from the BAM file: {}\n".format(args.extendReads,
                                                                                   int(frag_len_dict['median'])))

            fragment_length = args.fragmentLength

        else:
            # set as fragment length the read length
            fragment_length = int(read_len_dict['median'])
            if args.verbose:
                print "Estimated read length is {}".format(int(read_len_dict['median']))

        current_coverage = \
            float(bam_mapped * fragment_length) / args.normalizeTo1x
        # the scaling sets the coverage to match 1x
        scale_factor *= 1.0 / current_coverage
        if debug:
            print "Estimated current coverage {}".format(current_coverage)
            print "Scaling factor {}".format(args.scaleFactor)

    elif args.normalizeUsingRPKM:
        # the RPKM is the # reads per tile / \
        #    ( total reads (in millions) * tile length in Kb)
        millionReadsMapped = float(bam_mapped) / 1e6
        tileLengthInKb = float(args.binSize) / 1000

        scale_factor *= 1.0 / (millionReadsMapped * tileLengthInKb)

        if debug:
            print "scale factor using RPKM is {0}".format(args.scaleFactor)

    return scale_factor


def main(args=None):
    args = process_args(args)

    global debug
    if args.verbose:
        debug = 1
    else:
        debug = 0

    funcArgs = {'scaleFactor': get_scale_factor(args)}
    zerosToNans = True if args.missingDataAsZero == 'no' else False
    wr = writeBedGraph.WriteBedGraph([args.bam],
                                     binLength=args.binSize,
                                     stepSize=args.binSize,
                                     region=args.region,
                                     numberOfProcessors=args.numberOfProcessors,
                                     extendReads=args.extendReads,
                                     minMappingQuality=args.minMappingQuality,
                                     ignoreDuplicates=args.ignoreDuplicates,
                                     center_read=args.centerReads,
                                     zerosToNans=zerosToNans,
                                     samFlag_include=args.samFlagInclude,
                                     samFlag_exclude=args.samFlagExclude,
                                     verbose=args.verbose
                                     )

    wr.run(writeBedGraph.scaleCoverage, funcArgs, args.outFileName,
           format=args.outFileFormat, smooth_length=args.smoothLength)
