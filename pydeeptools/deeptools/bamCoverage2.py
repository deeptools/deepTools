import argparse
from deeptools import parserCommon
from deeptools.hp import r_bamcoverage


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
            usage='bamCoverage -b reads.bam -o coverage.bw\n'
            'help: bamCoverage -h / bamCoverage --help',
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
    if not args.smoothLength:
        args.smoothLength = 0
    if args.smoothLength and args.smoothLength <= args.binSize:
        print("Warning: the smooth length given ({}) is smaller than the bin "
              "size ({}).\n\n No smoothing will be done".format(args.smoothLength, args.binSize))
        args.smoothLength = 0

    if not args.ignoreForNormalization:
        args.ignoreForNormalization = []

    return args

def main(args=None):
    args = process_args(args)
    if not args.effectiveGenomeSize:
        args.effectiveGenomeSize = 0
    if not args.normalizeUsing:
        args.normalizeUsing = 'None'
    if not args.Offset:
        args.Offset = [1, -1]
    if not args.extendReads:
        args.extendReads = 0
    if not args.filterRNAstrand:
        args.filterRNAstrand = 'None'
    if not args.blackListFileName:
        args.blackListFileName = 'None'
    if not args.minMappingQuality:
        args.minMappingQuality = 0
    if not args.samFlagInclude:
        args.samFlagInclude = 0
    if not args.samFlagExclude:
        args.samFlagExclude = 0

    print(args)

    r_bamcoverage(
        args.bam, # bam file
        args.outFileName, # output file
        args.outFileFormat, # output format
        ## Normalization options
        args.normalizeUsing, # normalization
        args.effectiveGenomeSize, # effective genome size
        args.scaleFactor,
        # processing options
        args.MNase,
        args.Offset,
        args.extendReads,
        args.centerReads,
        args.filterRNAstrand,
        args.blackListFileName,
        args.ignoreForNormalization,
        args.skipNonCoveredRegions,
        args.smoothLength,
        args.binSize, # bin size
        # Filtering options
        args.ignoreDuplicates,
        args.minMappingQuality,
        args.samFlagInclude,
        args.samFlagExclude,
        args.minFragmentLength,
        args.maxFragmentLength,
        # running options
        args.numberOfProcessors, # threads
        [], # regions
        args.verbose # verbose
    )