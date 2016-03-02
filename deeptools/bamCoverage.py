#!/usr/bin/env python
# -*- coding: utf-8 -*-

# own tools
import argparse
from deeptools import writeBedGraph  # This should be made directly into a bigWig
from deeptools import parserCommon
from deeptools import bamHandler

debug = 0


def parseArguments():
    parentParser = parserCommon.getParentArgParse()
    bamParser = parserCommon.read_options()
    normalizationParser = parserCommon.normalization_options()
    requiredArgs = get_required_args()
    optionalArgs = get_optional_args()
    outputParser = parserCommon.output()
    parser = \
        argparse.ArgumentParser(
            parents=[requiredArgs, outputParser, optionalArgs,
                     parentParser, normalizationParser, bamParser],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            description='This tool takes an alignment of reads or fragments '
            'as input (BAM file) and generates a coverage track (bigWig or '
            'bedGraph) as output. '
            'The coverage is calculated as the number of reads per bin, '
            'where bins are short consecutive counting windows of a defined '
            'size. It is possible to extended the length of the reads '
            'to better reflect the actual fragment length. *bamCoverage* '
            'offers normalization by scaling factor, Reads Per Kilobase per '
            'Million mapped reads (RPKM), and 1x depth (reads per genome '
            'coverage, RPGC).\n',
            usage='An example usage is:'
            '$ bamCoverage -b reads.bam -o coverage.bw',
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

    optional.add_argument('--scaleFactor',
                          help='Indicate a number that you would like to use. When used in combination '
                          'with --normalizeTo1x or --normalizeUsingRPKM, the computed '
                          'scaling factor will be multiplied by the given scale factor.',
                          default=1.0,
                          type=float,
                          required=False)

    optional.add_argument('--MNase',
                          help='Determine nucleosome positions from MNase-seq data. '
                          'Only 3 nucleotides at the center of each fragment are counted. '
                          'The fragment ends are defined by the two mate reads. Only fragment lengths'
                          'between 130 - 200 bp are considered to avoid dinucleosomes or other artifacts.'
                          '*NOTE*: Requires paired-end data. A bin size of 1 is recommended.',
                          action='store_true')

    optional.add_argument('--filterRNAstrand',
                          help='Selects RNA-seq reads (single-end or paired-end) in '
                               'the given strand.',
                          choices=['forward', 'reverse'],
                          default=None)

    return parser


def scaleFactor(string):
    try:
        scalefactor1, scalefactor2 = string.split(":")
        scalefactors = (float(scalefactor1), float(scalefactor2))
    except:
        raise argparse.ArgumentTypeError(
            "Format of scaleFactors is factor1:factor2. "
            "The value given ( {} ) is not valid".format(string))

    return scalefactors


def process_args(args=None):
    args = parseArguments().parse_args(args)

    if args.scaleFactor != 1:
        args.normalizeTo1x = None
    if args.smoothLength and args.smoothLength <= args.binSize:
        print "Warning: the smooth length given ({}) is smaller than the bin "\
            "size ({}).\n\n No smoothing will be done".format(args.smoothLength, args.binSize)
        args.smoothLength = None

    return args


def get_scale_factor(args):

    scale_factor = args.scaleFactor
    bam_handle = bamHandler.openBam(args.bam)
    bam_mapped = parserCommon.bam_total_reads(bam_handle, args.ignoreForNormalization)
    blacklisted = parserCommon.bam_blacklisted_reads(bam_handle, args.ignoreForNormalization, args.blackListFileName)
    bam_mapped -= blacklisted

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
                    print("Fragment length based on paired en data "
                          "estimated to be {}".format(frag_len_dict['median']))

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
        million_reads_mapped = float(bam_mapped) / 1e6
        tile_len_in_kb = float(args.binSize) / 1000

        scale_factor *= 1.0 / (million_reads_mapped * tile_len_in_kb)

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

    func_args = {'scaleFactor': get_scale_factor(args)}

    if args.MNase:
        # check that library is paired end
        # using getFragmentAndReadSize
        from deeptools.getFragmentAndReadSize import get_read_and_fragment_length
        frag_len_dict, read_len_dict = get_read_and_fragment_length(args.bam,
                                                                    return_lengths=False,
                                                                    blackListFileName=args.blackListFileName,
                                                                    numberOfProcessors=args.numberOfProcessors,
                                                                    verbose=args.verbose)
        if frag_len_dict is None:
            exit("*Error*: For the --MNAse function a paired end library is required. ")

        wr = CenterFragment([args.bam],
                            binLength=args.binSize,
                            stepSize=args.binSize,
                            region=args.region,
                            blackListFileName=args.blackListFileName,
                            numberOfProcessors=args.numberOfProcessors,
                            extendReads=args.extendReads,
                            minMappingQuality=args.minMappingQuality,
                            ignoreDuplicates=args.ignoreDuplicates,
                            center_read=args.centerReads,
                            zerosToNans=args.skipNonCoveredRegions,
                            samFlag_include=args.samFlagInclude,
                            samFlag_exclude=args.samFlagExclude,
                            verbose=args.verbose,
                            )

    elif args.filterRNAstrand:
        wr = filterRnaStrand([args.bam],
                             binLength=args.binSize,
                             stepSize=args.binSize,
                             region=args.region,
                             numberOfProcessors=args.numberOfProcessors,
                             extendReads=args.extendReads,
                             minMappingQuality=args.minMappingQuality,
                             ignoreDuplicates=args.ignoreDuplicates,
                             center_read=args.centerReads,
                             zerosToNans=args.skipNonCoveredRegions,
                             samFlag_include=args.samFlagInclude,
                             samFlag_exclude=args.samFlagExclude,
                             verbose=args.verbose,
                             )

        wr.filter_strand = args.filterRNAstrand
    else:
        wr = writeBedGraph.WriteBedGraph([args.bam],
                                         binLength=args.binSize,
                                         stepSize=args.binSize,
                                         region=args.region,
                                         blackListFileName=args.blackListFileName,
                                         numberOfProcessors=args.numberOfProcessors,
                                         extendReads=args.extendReads,
                                         minMappingQuality=args.minMappingQuality,
                                         ignoreDuplicates=args.ignoreDuplicates,
                                         center_read=args.centerReads,
                                         zerosToNans=args.skipNonCoveredRegions,
                                         samFlag_include=args.samFlagInclude,
                                         samFlag_exclude=args.samFlagExclude,
                                         verbose=args.verbose,
                                         )

    wr.run(writeBedGraph.scaleCoverage, func_args, args.outFileName,
           blackListFileName=args.blackListFileName,
           format=args.outFileFormat, smoothLength=args.smoothLength)


class CenterFragment(writeBedGraph.WriteBedGraph):
    """
    Class to redefine the get_fragment_from_read for the --MNase case

    The coverage of the fragment is defined as the 2 or 3 basepairs at the
    center of the fragment length. d
    """
    def get_fragment_from_read(self, read):
        """
        Takes a proper pair fragment of high quality and limited
        to a certain length and outputs the center
        """
        fragment_start = fragment_end = None

        # only paired forward reads are considered
        # that are about one nuclesome turn long (130 - 200 bp)
        # TODO: this size range could be an input parameter.
        if read.is_proper_pair and not read.is_reverse and 130 < abs(read.tlen) < 200:
            # distance between pairs is even return two bases at the center
            if read.tlen % 2 == 0:
                fragment_start = read.pos + read.tlen / 2 - 1
                fragment_end = fragment_start + 2

            # distance is odd return three bases at the center
            else:
                fragment_start = read.pos + read.tlen / 2 - 1
                fragment_end = fragment_start + 3

        return [(fragment_start, fragment_end)]


class filterRnaStrand(writeBedGraph.WriteBedGraph):
    """
    Class to redefine the get_fragment_from_read for the --filterRNAstrand case

    Only reads either forward or reverse are kept as follows:

    For paired-end
    --------------
    reads forward:

     1. alignments of the second in pair (128) if they map to the forward strand (~16)
     2. alignments of the first in pair (64) if they map to the reverse  strand (~32)

     1. include 128, exclude 16
     or
     2. include 64 exclude 32

    reads reverse:
    1. alignments of the second in pair (128) if it maps to the reverse strand (16) 128 & 16 = 144
    2. alignments of the first in pair (64) if their mates map to the reverse strand (32) 64 & 32 = 96

     1. include 144
     or
     2. include 96

    For single-end
    --------------
    forward: include 16 (map forward strand)
    reverse: exclude 16

    """

    def get_fragment_from_read(self, read):
        """
        Gets only reads for the given strand
        """
        fragment_start = fragment_end = None

        # only paired forward reads are considered
        if read.is_paired:
            if self.filter_strand == 'forward':
                if (read.flag & 128 == 128 and read.flag & 16 == 0) or (read.flag & 64 == 64 and read.flag & 32 == 0):
                    return read.get_blocks()
            else:
                if read.flag & 144 == 144 or read.flag & 96 == 96:
                    return read.get_blocks()
        else:
            if self.filter_strand == 'forward' and read.flag & 16 == 16:
                return read.get_blocks()
            elif self.filter_strand == 'reverse' and read.flag & 16 == 0:
                return read.get_blocks()

        return [(fragment_start, fragment_end)]
