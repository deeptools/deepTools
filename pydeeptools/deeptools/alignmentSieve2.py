import argparse
import signal
import sys
from importlib.metadata import version

from deeptools import parserCommon
from deeptools.hp import r_alignmentsieve


def parseArguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="This tool filters alignments in a BAM/CRAM file according the the specified parameters. It can optionally output to BEDPE format.",
        usage='alignmentSieve -b sample1.bam -o sample1.filtered.bam --minMappingQuality 10 --filterMetrics log.txt\n'
        'help: alignmentSieve -h / alignmentSieve --help')

    required = parser.add_argument_group('Required arguments')
    required.add_argument('--bam', '-b',
                          metavar='FILE1',
                          help='An indexed BAM file.',
                          required=True)

    required.add_argument('--outFile', '-o',
                          help='The file to write results to. These are the alignments or fragments that pass the filtering criteria.')

    general = parser.add_argument_group('General arguments')
    general.add_argument('--numberOfProcessors', '-p',
                         help='Number of processors to use. Type "max/2" to '
                         'use half the maximum number of processors or "max" '
                         'to use all available processors. (Default: %(default)s)',
                         metavar="INT",
                         type=parserCommon.numberOfProcessors,
                         default=1,
                         required=False)

    general.add_argument('--filterMetrics',
                         metavar="FILE.log",
                         default="None",
                         help="The number of entries in total and filtered are saved to this file")

    general.add_argument('--filteredOutReads',
                          metavar="filtered.bam",
                          default="None",
                          help="If desired, all reads NOT passing the filtering criteria can be written to this file.")

    general.add_argument('--label', '-l',
                          metavar='sample1',
                          default="None",
                          help='User defined label instead of the default label '
                          '(file name).')

    general.add_argument('--smartLabels',
                          action='store_true',
                          help='Instead of manually specifying a labels for the input '
                          'file, this causes deepTools to use the file name '
                          'after removing the path and extension.')

    general.add_argument('--verbose', '-v',
                         help='Set to see processing messages.',
                         action='store_true')

    general.add_argument('--version', action='version',
                         version='%(prog)s {}'.format(version('deeptools')))

    general.add_argument('--shift',
                         nargs='+',
                         type=int,
                         help='Shift the left and right end of a read (for BAM files) or a fragment (for BED files). A positive value shift an end to the right (on the + strand) and a negative value shifts a fragment to the left. Either 2 or 4 integers can be provided. For example, "2 -3" will shift the left-most fragment end two bases to the right and the right-most end 3 bases to the left. If 4 integers are provided, then the first and last two refer to fragments whose read 1 is on the left or right, respectively. Consequently, it is possible to take strand into consideration for strand-specific protocols. A fragment whose length falls below 1 due to shifting will not be written to the output. See the online documentation for graphical examples. Note that non-properly-paired reads will be filtered.')

    general.add_argument('--ATACshift',
                         action='store_true',
                         help='Shift the produced BAM file or BEDPE regions as commonly done for ATAC-seq. This is equivalent to --shift 4 -5 5 -4.')

    output = parser.add_argument_group('Output arguments')
    output.add_argument('--BED',
                        action='store_true',
                        help='Instead of producing BAM files, write output in BEDPE format (as defined by MACS2). Note that only reads/fragments passing filtering criterion are written in BEDPE format.')

    filtering = parser.add_argument_group('Optional arguments')

    filtering.add_argument('--filterRNAstrand',
                           help='Selects RNA-seq reads (single-end or paired-end) in '
                                'the given strand. (Default: %(default)s)',
                           choices=['forward', 'reverse', 'None'],
                           default='None')

    filtering.add_argument('--minMappingQuality',
                           metavar='INT',
                           help='If set, only reads that have a mapping '
                           'quality score of at least this are '
                           'considered.',
                           default=0,
                           type=int)

    filtering.add_argument('--samFlagInclude',
                           help='Include reads based on the SAM flag. For example, '
                           'to get only reads that are the first mate, use a flag of 64. '
                           'This is useful to count properly paired reads only once, '
                           'as otherwise the second mate will be also considered for the '
                           'coverage.',
                           metavar='INT',
                           type=int,
                           default=0,
                           required=False)

    filtering.add_argument('--samFlagExclude',
                           help='Exclude reads based on the SAM flag. For example, '
                           'to get only reads that map to the forward strand, use '
                           '--samFlagExclude 16, where 16 is the SAM flag for reads '
                           'that map to the reverse strand.',
                           metavar='INT',
                           default=0,
                           type=int,
                           required=False)

    filtering.add_argument('--blackListFileName', '-bl',
                           help="A BED or GTF file (optionally gzip-compressed) containing regions that should be excluded from all analyses. Filtering is performed at base-pair resolution, so only the portion of a read/fragment that overlaps a blacklisted region is excluded. Please note that you should adjust the effective genome size, if relevant.",
                           metavar="BED file",
                           nargs="+",
                           default="None",
                           required=False)

    filtering.add_argument('--ignoreDuplicates',
                            help='If set, reads that are marked as PCR '
                            'or optical duplicates (SAM flag 0x400) will '
                            'be filtered out.',
                            action='store_true')

    filtering.add_argument('--minFragmentLength',
                           help='The minimum fragment length needed for read/pair '
                           'inclusion. This option is primarily useful '
                           'in ATACseq experiments, for filtering mono- or '
                           'di-nucleosome fragments. (Default: %(default)s)',
                           metavar='INT',
                           default=0,
                           type=int,
                           required=False)

    filtering.add_argument('--maxFragmentLength',
                           help='The maximum fragment length needed for read/pair '
                           'inclusion. A value of 0 indicates no limit. (Default: %(default)s)',
                           metavar='INT',
                           default=0,
                           type=int,
                           required=False)

    return parser


def main(args=None):
    args = parseArguments().parse_args(args)
    if args.shift:
        if len(args.shift) not in [2, 4]:
            sys.exit("The --shift option can accept either 2 or 4 values only.")
        if len(args.shift) == 2:
            args.shift.extend([-args.shift[1], -args.shift[0]])
    else:
        args.shift = []
    if args.ATACshift:
        if args.shift:
            print("Warning! The --ATACshift option is used, but a --shift option is provided as well. The latter will be ignored in favor of 4 -5 5 -4.")
        args.shift = [4, -5, 5, -4]

    if args.ignoreDuplicates:
        args.samFlagExclude |= 0x400

    if not args.blackListFileName:
        args.blackListFileName = "None"
    elif isinstance(args.blackListFileName, list):
        if len(args.blackListFileName) != 1:
            sys.exit("Only one blacklist file is supported when using '--alignmentsieve rust'.")
        args.blackListFileName = args.blackListFileName[0]

    signal.signal(signal.SIGINT, signal.SIG_DFL)
    r_alignmentsieve(
        args.bam,
        args.outFile,
        args.numberOfProcessors,
        args.filterMetrics,
        args.filteredOutReads,
        args.verbose,
        args.shift,
        args.BED,
        args.filterRNAstrand,
        args.minMappingQuality,
        args.samFlagInclude,
        args.samFlagExclude,
        args.blackListFileName,
        args.minFragmentLength,
        args.maxFragmentLength,
        args.label,
        args.smartLabels,
    )
