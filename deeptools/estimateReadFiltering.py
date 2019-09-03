#!/usr/bin/env python
import argparse
import sys

from deeptools import parserCommon, bamHandler, utilities
from deeptools.mapReduce import mapReduce
from deeptools.utilities import smartLabels
from deeptools._version import __version__


def parseArguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
This tool estimates the number of reads that would be filtered given a set of
settings and prints this to the terminal. Further, it tracks the number of singleton reads. The following metrics will always be tracked regardless of what you specify (the order output also matches this):

 * Total reads (including unmapped)
 * Mapped reads
 * Reads in blacklisted regions (--blackListFileName)

The following metrics are estimated according to the --binSize and --distanceBetweenBins parameters
 * Estimated mapped reads filtered (the total number of mapped reads filtered for any reason)
 * Alignments with a below threshold MAPQ (--minMappingQuality)
 * Alignments with at least one missing flag (--samFlagInclude)
 * Alignments with undesirable flags (--samFlagExclude)
 * Duplicates determined by deepTools (--ignoreDuplicates)
 * Duplicates marked externally (e.g., by picard)
 * Singletons (paired-end reads with only one mate aligning)
 * Wrong strand (due to --filterRNAstrand)

The sum of these may be more than the total number of reads. Note that alignments are sampled from bins of size --binSize spaced --distanceBetweenBins apart.
""",
        usage='Example usage: estimateReadFiltering.py -b sample1.bam sample2.bam > log.txt')

    required = parser.add_argument_group('Required arguments')
    required.add_argument('--bamfiles', '-b',
                          metavar='FILE1 FILE2',
                          help='List of indexed bam files separated by spaces.',
                          nargs='+',
                          required=True)

    general = parser.add_argument_group('General arguments')

    general.add_argument('--outFile', '-o',
                         type=parserCommon.writableFile,
                         help='The file to write results to. By default, results are printed to the console')

    general.add_argument('--sampleLabels',
                         help='Labels for the samples. The '
                         'default is to use the file name of the '
                         'sample. The sample labels should be separated '
                         'by spaces and quoted if a label itself'
                         'contains a space E.g. --sampleLabels label-1 "label 2"  ',
                         nargs='+')

    general.add_argument('--smartLabels',
                         action='store_true',
                         help='Instead of manually specifying labels for the input '
                         'BAM files, this causes deepTools to use the '
                         'file name after removing the path and extension.')

    general.add_argument('--binSize', '-bs',
                         metavar='INT',
                         help='Length in bases of the window used to sample the genome. (Default: %(default)s)',
                         default=1000000,
                         type=int)

    general.add_argument('--distanceBetweenBins', '-n',
                         metavar='INT',
                         help='To reduce the computation time, not every possible genomic '
                         'bin is sampled. This option allows you to set the distance '
                         'between bins actually sampled from. Larger numbers are sufficient '
                         'for high coverage samples, while smaller values are useful for '
                         'lower coverage samples. Note that if you specify a value that '
                         'results in too few (<1000) reads sampled, the value will be '
                         'decreased. (Default: %(default)s)',
                         default=10000,
                         type=int)

    general.add_argument('--numberOfProcessors', '-p',
                         help='Number of processors to use. Type "max/2" to '
                         'use half the maximum number of processors or "max" '
                         'to use all available processors. (Default: %(default)s)',
                         metavar="INT",
                         type=parserCommon.numberOfProcessors,
                         default=1,
                         required=False)

    general.add_argument('--verbose', '-v',
                         help='Set to see processing messages.',
                         action='store_true')

    general.add_argument('--version', action='version',
                         version='%(prog)s {}'.format(__version__))

    filtering = parser.add_argument_group('Optional arguments')

    filtering.add_argument('--filterRNAstrand',
                           help='Selects RNA-seq reads (single-end or paired-end) in '
                                'the given strand. (Default: %(default)s)',
                           choices=['forward', 'reverse'],
                           default=None)

    filtering.add_argument('--ignoreDuplicates',
                           help='If set, reads that have the same orientation '
                           'and start position will be considered only '
                           'once. If reads are paired, the mate\'s position '
                           'also has to coincide to ignore a read.',
                           action='store_true')

    filtering.add_argument('--minMappingQuality',
                           metavar='INT',
                           help='If set, only reads that have a mapping '
                           'quality score of at least this are '
                           'considered.',
                           type=int)

    filtering.add_argument('--samFlagInclude',
                           help='Include reads based on the SAM flag. For example, '
                           'to get only reads that are the first mate, use a flag of 64. '
                           'This is useful to count properly paired reads only once, '
                           'as otherwise the second mate will be also considered for the '
                           'coverage. (Default: %(default)s)',
                           metavar='INT',
                           default=None,
                           type=int,
                           required=False)

    filtering.add_argument('--samFlagExclude',
                           help='Exclude reads based on the SAM flag. For example, '
                           'to get only reads that map to the forward strand, use '
                           '--samFlagExclude 16, where 16 is the SAM flag for reads '
                           'that map to the reverse strand. (Default: %(default)s)',
                           metavar='INT',
                           default=None,
                           type=int,
                           required=False)

    filtering.add_argument('--blackListFileName', '-bl',
                           help="A BED or GTF file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant.",
                           metavar="BED file",
                           nargs="+",
                           required=False)

    return parser


def getFiltered_worker(arglist):
    chrom, start, end, args = arglist
    # Fix the bounds
    if end - start > args.binSize and end - start > args.distanceBetweenBins:
        end -= args.distanceBetweenBins
    if end <= start:
        end = start + 1

    o = []
    for fname in args.bamfiles:
        fh = bamHandler.openBam(fname)
        chromUse = utilities.mungeChromosome(chrom, fh.references)
        prev_pos = set()
        lpos = None

        minMapq = 0
        samFlagInclude = 0
        samFlagExclude = 0
        internalDupes = 0
        externalDupes = 0
        singletons = 0
        filterRNAstrand = 0
        nFiltered = 0
        total = 0  # This is only used to estimate the percentage affected
        for read in fh.fetch(chromUse, start, end):
            filtered = 0
            if read.pos < start:
                # ensure that we never double count (in case distanceBetweenBins == 0)
                continue

            if read.flag & 4:
                # Ignore unmapped reads, they were counted already
                continue

            if args.minMappingQuality and read.mapq < args.minMappingQuality:
                filtered = 1
                minMapq += 1
            if args.samFlagInclude and read.flag & args.samFlagInclude != args.samFlagInclude:
                filtered = 1
                samFlagInclude += 1
            if args.samFlagExclude and read.flag & args.samFlagExclude != 0:
                filtered = 1
                samFlagExclude += 1
            if args.ignoreDuplicates:
                # Assuming more or less concordant reads, use the fragment bounds, otherwise the start positions
                if read.tlen >= 0:
                    s = read.pos
                    e = s + read.tlen
                else:
                    s = read.pnext
                    e = s - read.tlen
                if read.reference_id != read.next_reference_id:
                    e = read.pnext
                if lpos is not None and lpos == read.reference_start \
                        and (s, e, read.next_reference_id, read.is_reverse) in prev_pos:
                    filtered = 1
                    internalDupes += 1
                if lpos != read.reference_start:
                    prev_pos.clear()
                lpos = read.reference_start
                prev_pos.add((s, e, read.next_reference_id, read.is_reverse))
            if read.is_duplicate:
                filtered = 1
                externalDupes += 1
            if read.is_paired and read.mate_is_unmapped:
                filtered = 1
                singletons += 1

            # filterRNAstrand
            if args.filterRNAstrand:
                if read.is_paired:
                    if args.filterRNAstrand == 'forward':
                        if read.flag & 144 == 128 or read.flag & 96 == 64:
                            pass
                        else:
                            filtered = 1
                            filterRNAstrand += 1
                    elif args.filterRNAstrand == 'reverse':
                        if read.flag & 144 == 144 or read.flag & 96 == 96:
                            pass
                        else:
                            filtered = 1
                            filterRNAstrand += 1
                else:
                    if args.filterRNAstrand == 'forward':
                        if read.flag & 16 == 16:
                            pass
                        else:
                            filtered = 1
                            filterRNAstrand += 1
                    elif args.filterRNAstrand == 'reverse':
                        if read.flag & 16 == 0:
                            pass
                        else:
                            filtered = 1
                            filterRNAstrand += 1

            total += 1
            nFiltered += filtered
        fh.close()

        # Append a tuple to the output
        tup = (total, nFiltered, minMapq, samFlagInclude, samFlagExclude, internalDupes, externalDupes, singletons, filterRNAstrand)
        o.append(tup)
    return o


def main(args=None):
    args = parseArguments().parse_args(args)

    if not args.sampleLabels and args.smartLabels:
        args.sampleLabels = smartLabels(args.bamfiles)

    if args.sampleLabels and len(args.sampleLabels) != len(args.bamfiles):
        sys.stderr.write("\nError: --sampleLabels specified but it doesn't match the number of BAM files!\n")
        sys.exit(1)

    if args.outFile is None:
        of = sys.stdout
    else:
        of = open(args.outFile, "w")

    bhs = [bamHandler.openBam(x, returnStats=True, nThreads=args.numberOfProcessors) for x in args.bamfiles]
    mapped = [x[1] for x in bhs]
    unmappedList = [x[2] for x in bhs]
    bhs = [x[0] for x in bhs]

    # Get the reads in blacklisted regions
    if args.blackListFileName:
        blacklisted = []
        for bh in bhs:
            blacklisted.append(utilities.bam_blacklisted_reads(bh, None, args.blackListFileName, args.numberOfProcessors))
    else:
        blacklisted = [0] * len(bhs)

    # Get the total and mapped reads
    total = [x + y for x, y in list(zip(mapped, unmappedList))]

    chrom_sizes = list(zip(bhs[0].references, bhs[0].lengths))
    for x in bhs:
        x.close()

    # Get the remaining metrics
    res = mapReduce([args],
                    getFiltered_worker,
                    chrom_sizes,
                    genomeChunkLength=args.binSize + args.distanceBetweenBins,
                    blackListFileName=args.blackListFileName,
                    numberOfProcessors=args.numberOfProcessors,
                    verbose=args.verbose)

    totals = [0] * len(args.bamfiles)
    nFiltered = [0] * len(args.bamfiles)
    MAPQs = [0] * len(args.bamfiles)
    flagIncludes = [0] * len(args.bamfiles)
    flagExcludes = [0] * len(args.bamfiles)
    internalDupes = [0] * len(args.bamfiles)
    externalDupes = [0] * len(args.bamfiles)
    singletons = [0] * len(args.bamfiles)
    rnaStrand = [0] * len(args.bamfiles)
    for x in res:
        for idx, r in enumerate(x):
            totals[idx] += r[0]
            nFiltered[idx] += r[1]
            MAPQs[idx] += r[2]
            flagIncludes[idx] += r[3]
            flagExcludes[idx] += r[4]
            internalDupes[idx] += r[5]
            externalDupes[idx] += r[6]
            singletons[idx] += r[7]
            rnaStrand[idx] += r[8]

    # Print some output
    of.write("Sample\tTotal Reads\tMapped Reads\tAlignments in blacklisted regions\tEstimated mapped reads filtered\tBelow MAPQ\tMissing Flags\tExcluded Flags\tInternally-determined Duplicates\tMarked Duplicates\tSingletons\tWrong strand\n")
    for idx, _ in enumerate(args.bamfiles):
        if args.sampleLabels:
            of.write(args.sampleLabels[idx])
        else:
            of.write(args.bamfiles[idx])
        of.write("\t{}\t{}\t{}".format(total[idx], mapped[idx], blacklisted[idx]))
        # nFiltered
        metric = 0.0
        if totals[idx] > 0:
            metric = blacklisted[idx] + float(nFiltered[idx]) / float(totals[idx]) * mapped[idx]
        of.write("\t{}".format(min(round(metric, 1), mapped[idx])))
        # MAPQ
        metric = 0.0
        if totals[idx] > 0:
            metric = float(MAPQs[idx]) / float(totals[idx]) * mapped[idx]
        of.write("\t{}".format(min(round(metric, 1), mapped[idx])))
        # samFlagInclude
        metric = 0.0
        if totals[idx] > 0:
            metric = float(flagIncludes[idx]) / float(totals[idx]) * mapped[idx]
        of.write("\t{}".format(min(round(metric, 1), mapped[idx])))
        # samFlagExclude
        metric = 0.0
        if totals[idx] > 0:
            metric = float(flagExcludes[idx]) / float(totals[idx]) * mapped[idx]
        of.write("\t{}".format(min(round(metric, 1), mapped[idx])))
        # Internally determined duplicates
        metric = 0.0
        if totals[idx] > 0:
            metric = float(internalDupes[idx]) / float(totals[idx]) * mapped[idx]
        of.write("\t{}".format(min(round(metric, 1), mapped[idx])))
        # Externally marked duplicates
        metric = 0.0
        if totals[idx] > 0:
            metric = float(externalDupes[idx]) / float(totals[idx]) * mapped[idx]
        of.write("\t{}".format(min(round(metric, 1), mapped[idx])))
        # Singletons
        metric = 0.0
        if totals[idx] > 0:
            metric = float(singletons[idx]) / float(totals[idx]) * mapped[idx]
        of.write("\t{}".format(min(round(metric, 1), mapped[idx])))
        # filterRNAstrand
        metric = 0.0
        if totals[idx] > 0:
            metric = float(rnaStrand[idx]) / float(totals[idx]) * mapped[idx]
        of.write("\t{}".format(min(round(metric, 1), mapped[idx])))
        of.write("\n")

    if args.outFile is not None:
        of.close()

    return 0
