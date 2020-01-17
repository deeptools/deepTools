#!/usr/bin/env python
# -*- coding: utf-8 -*-

# own tools
import argparse
import sys
import numpy as np
from deeptools import writeBedGraph  # This should be made directly into a bigWig
from deeptools import parserCommon
from deeptools.getScaleFactor import get_scale_factor
from deeptools.bamHandler import openBam

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
            'Million mapped reads (RPKM), counts per million (CPM), bins per '
            'million mapped reads (BPM) and 1x depth (reads per genome '
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
                          help='The computed scaling factor (or 1, if not applicable) will '
                          'be multiplied by this. (Default: %(default)s)',
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
                          'offset of 0 is not permitted. If two values are specified, then they will be '
                          'used to specify a range of positions. Note that specifying something like '
                          '--Offset 5 -1 will result in the 5th through last position being used, which '
                          'is equivalent to trimming 4 bases from the 5-prime end of alignments. Note '
                          'that if you specify --centerReads, the centering will be performed before the '
                          'offset.',
                          metavar='INT',
                          type=int,
                          nargs='+',
                          required=False)

    optional.add_argument('--filterRNAstrand',
                          help='Selects RNA-seq reads (single-end or paired-end) originating from genes '
                          'on the given strand. This option assumes a standard dUTP-based library '
                          'preparation (that is, --filterRNAstrand=forward keeps minus-strand reads, '
                          'which originally came from genes on the forward strand using a dUTP-based '
                          'method). Consider using --samExcludeFlag instead for filtering by strand in '
                          'other contexts.',
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

    if args.smoothLength and args.smoothLength <= args.binSize:
        print("Warning: the smooth length given ({}) is smaller than the bin "
              "size ({}).\n\n No smoothing will be done".format(args.smoothLength, args.binSize))
        args.smoothLength = None

    if not args.ignoreForNormalization:
        args.ignoreForNormalization = []

    return args


def main(args=None):
    args = process_args(args)

    global debug
    if args.verbose:
        sys.stderr.write("Specified --scaleFactor: {}\n".format(args.scaleFactor))
        debug = 1
    else:
        debug = 0

    if args.normalizeUsing == 'None':
        args.normalizeUsing = None  # For the sake of sanity
    elif args.normalizeUsing == 'RPGC' and not args.effectiveGenomeSize:
        sys.exit("RPGC normalization requires an --effectiveGenomeSize!\n")

    if args.normalizeUsing:
        # if a normalization is required then compute the scale factors
        bam, mapped, unmapped, stats = openBam(args.bam, returnStats=True, nThreads=args.numberOfProcessors)
        bam.close()
        scale_factor = get_scale_factor(args, stats)
    else:
        scale_factor = args.scaleFactor

    func_args = {'scaleFactor': scale_factor}

    # This fixes issue #520, where --extendReads wasn't honored if --filterRNAstrand was used
    if args.filterRNAstrand and not args.Offset:
        args.Offset = [1, -1]

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
            sys.exit("*Error*: For the --MNAse function a paired end library is required. ")

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
                            chrsToSkip=args.ignoreForNormalization,
                            verbose=args.verbose,
                            )

    elif args.Offset:
        if len(args.Offset) > 1:
            if args.Offset[0] == 0:
                sys.exit("*Error*: An offset of 0 isn't allowed, since offsets are 1-based positions inside each alignment.")
            if args.Offset[1] > 0 and args.Offset[1] < args.Offset[0]:
                sys.exir("'Error*: The right side bound is less than the left-side bound. This is inappropriate.")
        else:
            if args.Offset[0] == 0:
                sys.exit("*Error*: An offset of 0 isn't allowed, since offsets are 1-based positions inside each alignment.")
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
                            chrsToSkip=args.ignoreForNormalization,
                            verbose=args.verbose)
        wr.filter_strand = args.filterRNAstrand
        wr.Offset = args.Offset
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
                                         chrsToSkip=args.ignoreForNormalization,
                                         verbose=args.verbose,
                                         )

    wr.run(writeBedGraph.scaleCoverage, func_args, args.outFileName,
           blackListFileName=args.blackListFileName,
           format=args.outFileFormat, smoothLength=args.smoothLength)


class OffsetFragment(writeBedGraph.WriteBedGraph):
    """
    Class to redefine the get_fragment_from_read for the --Offset case
    """
    def filterStrand(self, read, rv):
        """
        A generic read filtering function that gets used by everything in this class.

        rv is returned if the strand is correct, otherwise [(None, None)]
        """
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

    def get_fragment_from_read_list(self, read, offset):
        """
        Return the range of exons from the 0th through 1st bases, inclusive. Positions are 1-based
        """
        rv = [(None, None)]
        blocks = read.get_blocks()
        blockLen = sum([x[1] - x[0] for x in blocks])

        if self.defaultFragmentLength != 'read length':
            if self.is_proper_pair(read, self.maxPairedFragmentLength):
                if read.is_reverse:
                    foo = (read.next_reference_start, read.reference_start)
                    if foo[0] < foo[1]:
                        blocks.insert(0, foo)
                else:
                    foo = (read.reference_end, read.reference_end + abs(read.template_length) - read.infer_query_length())
                    if foo[0] < foo[1]:
                        blocks.append(foo)

            # Extend using the default fragment length
            else:
                if read.is_reverse:
                    foo = (read.reference_start - self.defaultFragmentLength + read.infer_query_length(), read.reference_start)
                    if foo[0] < 0:
                        foo = (0, foo[1])
                    if foo[0] < foo[1]:
                        blocks.insert(0, foo)
                else:
                    foo = (read.reference_end, read.reference_end + self.defaultFragmentLength - read.infer_query_length())
                    if foo[0] < foo[1]:
                        blocks.append(foo)

        stretch = []
        # For the sake of simplicity, convert [(10, 20), (30, 40)] to [10, 11, 12, 13, ..., 40]
        # Then subset accordingly
        for block in blocks:
            stretch.extend(range(block[0], block[1]))
        if read.is_reverse:
            stretch = stretch[::-1]

        # Handle --centerReads
        if self.center_read:
            _ = (len(stretch) - blockLen) // 2
            stretch = stretch[_:_ + blockLen]

        # Subset by --Offset
        try:
            foo = stretch[offset[0]:offset[1]]
        except:
            return rv

        if len(foo) == 0:
            return rv
        if read.is_reverse:
            foo = foo[::-1]

        # Convert the stretch back to a list of tuples
        foo = np.array(foo)
        d = foo[1:] - foo[:-1]
        idx = np.argwhere(d > 1).flatten().tolist()  # This now holds the interval bounds as a list
        idx.append(-1)
        last = 0
        rv = []
        for i in idx:
            rv.append((foo[last].astype("int"), foo[i].astype("int") + 1))
            last = i + 1

        # Handle strand filtering, if needed
        return self.filterStrand(read, rv)

    def get_fragment_from_read(self, read):
        """
        This is mostly a wrapper for self.get_fragment_from_read_list(),
        which needs a list and for the offsets to be tweaked by 1.
        """
        offset = [x for x in self.Offset]
        if len(offset) > 1:
            if offset[0] > 0:
                offset[0] -= 1
            if offset[1] < 0:
                offset[1] += 1
        else:
            if offset[0] > 0:
                offset[0] -= 1
                offset = [offset[0], offset[0] + 1]
            else:
                if offset[0] < -1:
                    offset = [offset[0], offset[0] + 1]
                else:
                    offset = [offset[0], None]
        if offset[1] == 0:
            # -1 gets switched to 0, which screws things up
            offset = (offset[0], None)
        return self.get_fragment_from_read_list(read, offset)


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
