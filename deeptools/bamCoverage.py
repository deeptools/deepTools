#!/usr/bin/env python
# -*- coding: utf-8 -*-

# own tools
import argparse
from deeptools import writeBedGraph  # This should be made directly into a bigWig
from deeptools import parserCommon
from deeptools.getScaleFactor import get_scale_factor

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
                          'between 130 - 200 bp are considered to avoid dinucleosomes or other artifacts. '
                          'By default, any fragments smaller or larger than this are ignored. To '
                          'over-ride this, use the --minFragmentLength and --maxFragmentLength options, '
                          'which will default to 130 and 200 if not otherwise specified in the presence '
                          'of --MNase. *NOTE*: Requires paired-end data. A bin size of 1 is recommended.',
                          action='store_true')

    optional.add_argument('--Offset',
                          help='Uses this offset inside of each read as the signal. This is useful in '
                          'cases like RiboSeq or GROseq, where the signal is 12, 15 or 0 bases past the '
                          'start of the read. This can be paired with the --filterRNAstrand option. '
                          'Note that negative values indicate offsets from the end of each read. A value '
                          'of 1 indicates the first base of the alignment (taking alignment orientation '
                          'into account). Likewise, a value of -1 is the last base of the alignment. An '
                          'offset of 0 is not permitted.',
                          metavar='INT',
                          type=int,
                          required=False)

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
        print("Warning: the smooth length given ({}) is smaller than the bin "
              "size ({}).\n\n No smoothing will be done".format(args.smoothLength, args.binSize))
        args.smoothLength = None

    return args


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

        # Set some default fragment length bounds
        if args.minFragmentLength == 0:
            args.minFragmentLength = 130
        if args.maxFragmentLength == 0:
            args.maxFragmentLength = 200

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
                            minFragmentLength=args.minFragmentLength,
                            maxFragmentLength=args.maxFragmentLength,
                            verbose=args.verbose,
                            )

    elif args.Offset:
        if args.Offset == 0:
            exit("*Error*: An offset of 0 isn't allowed, since offsets are 1-based positions inside each alignment.")
        wr = OffsetFragment([args.bam],
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
                            minFragmentLength=args.minFragmentLength,
                            maxFragmentLength=args.maxFragmentLength,
                            verbose=args.verbose)
        wr.filter_strand = args.filterRNAstrand
        wr.Offset = args.Offset

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
                             minFragmentLength=args.minFragmentLength,
                             maxFragmentLength=args.maxFragmentLength,
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
                                         minFragmentLength=args.minFragmentLength,
                                         maxFragmentLength=args.maxFragmentLength,
                                         verbose=args.verbose,
                                         )

    wr.run(writeBedGraph.scaleCoverage, func_args, args.outFileName,
           blackListFileName=args.blackListFileName,
           format=args.outFileFormat, smoothLength=args.smoothLength)


class OffsetFragment(writeBedGraph.WriteBedGraph):
    """
    Class to redefine the get_fragment_from_read for the --Offset case
    """
    def get_fragment_from_read(self, read):
        rv = [(None, None)]
        if self.Offset > read.query_length:
            return rv
        if read.is_paired:
            return rv
        blocks = read.get_blocks()
        foo = self.Offset
        if foo < 0:
            foo = read.infer_query_length() + foo + 1
        if read.is_reverse:
            for idx in range(len(blocks)):
                block = blocks[-idx - 1]
                if block[1] - block[0] >= foo:
                    rv = [(block[1] - foo, block[1] - foo + 1)]
                    break
                foo -= block[1] - block[0]
        else:
            for block in blocks:
                if block[1] - block[0] >= foo:
                    rv = [(block[0] + foo - 1, block[0] + foo)]
                    break
                foo -= block[1] - block[0]

        # Filter by RNA strand, if desired
        if read.is_paired:
            if self.filter_strand == 'forward':
                if read.flag & 144 == 128 or read.flag & 96 == 64:
                    return rv
            elif self.filter_strand == 'reverse':
                if read.flag & 144 == 144 or read.flag & 96 == 96:
                    return rv
            else:
                return rv
        else:
            if self.filter_strand == 'forward':
                if read.flag & 16 == 16:
                    return rv
            elif self.filter_strand == 'reverse':
                if read.flag & 16 == 0:
                    return rv
            else:
                return rv

        return [(None, None)]


class CenterFragment(writeBedGraph.WriteBedGraph):
    """
    Class to redefine the get_fragment_from_read for the --MNase case

    The coverage of the fragment is defined as the 2 or 3 basepairs at the
    center of the fragment length.
    """
    def get_fragment_from_read(self, read):
        """
        Takes a proper pair fragment of high quality and limited
        to a certain length and outputs the center
        """
        fragment_start = fragment_end = None

        # only paired forward reads are considered
        # Fragments have already been filtered according to length
        if read.is_proper_pair and not read.is_reverse and 1 < abs(read.tlen):
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
