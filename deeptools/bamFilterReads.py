#!/usr/bin/env python
import argparse
import pysam
import os

from deeptools import parserCommon
from deeptools.bamHandler import openBam
from deeptools.mapReduce import mapReduce
from deeptools._version import __version__
from deeptools.utilities import getTLen, smartLabels, getTempFileName


def parseArguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="This tool filters alignments in a BAM/CRAM file according the the specified parameters.",
        usage='Example usage: bamFilterReads.py -b sample1.bam -o sample1.filtered.bam --minMappingQuality 10 --filterMetrics log.txt')

    required = parser.add_argument_group('Required arguments')
    required.add_argument('--bam', '-b',
                          metavar='FILE1',
                          help='An indexed BAM file.',
                          required=True)

    required.add_argument('--outFile', '-o',
                          help='The file to write results to.')

    general = parser.add_argument_group('General arguments')
    general.add_argument('--numberOfProcessors', '-p',
                         help='Number of processors to use. Type "max/2" to '
                         'use half the maximum number of processors or "max" '
                         'to use all available processors.',
                         metavar="INT",
                         type=parserCommon.numberOfProcessors,
                         default=1,
                         required=False)

    general.add_argument('--filterMetrics',
                         metavar="FILE.log",
                         help="The number of entries in total and filtered are saved to this file")

    general.add_argument('--label', '-l',
                         metavar='sample1',
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
                         version='%(prog)s {}'.format(__version__))

    filtering = parser.add_argument_group('Optional arguments')

    filtering.add_argument('--filterRNAstrand',
                           help='Selects RNA-seq reads (single-end or paired-end) in '
                                'the given strand.',
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
                           'coverage.',
                           metavar='INT',
                           default=None,
                           type=int,
                           required=False)

    filtering.add_argument('--samFlagExclude',
                           help='Exclude reads based on the SAM flag. For example, '
                           'to get only reads that map to the forward strand, use '
                           '--samFlagExclude 16, where 16 is the SAM flag for reads '
                           'that map to the reverse strand.',
                           metavar='INT',
                           default=None,
                           type=int,
                           required=False)

    filtering.add_argument('--blackListFileName', '-bl',
                           help="A BED or GTF file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant.",
                           metavar="BED file",
                           nargs="+",
                           required=False)

    filtering.add_argument('--minFragmentLength',
                           help='The minimum fragment length needed for read/pair '
                           'inclusion. This option is primarily useful '
                           'in ATACseq experiments, for filtering mono- or '
                           'di-nucleosome fragments.',
                           metavar='INT',
                           default=0,
                           type=int,
                           required=False)

    filtering.add_argument('--maxFragmentLength',
                           help='The maximum fragment length needed for read/pair '
                           'inclusion.',
                           metavar='INT',
                           default=0,
                           type=int,
                           required=False)

    return parser


def filterWorker(arglist):
    chrom, start, end, args = arglist
    fh = openBam(args.bam)

    oname = getTempFileName(suffix='bam')
    ofh = pysam.AlignmentFile(oname, mode='wbu', template=fh)

    prev_pos = set()
    lpos = None

    nFiltered = 0
    total = 0
    for read in fh.fetch(chrom, start, end):
        if read.pos < start:
            # ensure that we never double count (in case distanceBetweenBins == 0)
            continue

        total += 1
        if read.flag & 4:
            # Ignore unmapped reads, they were counted already
            nFiltered += 1
            continue

        if args.minMappingQuality and read.mapq < args.minMappingQuality:
            nFiltered += 1
            continue

        if args.samFlagInclude and read.flag & args.samFlagInclude != args.samFlagInclude:
            nFiltered += 1
            continue

        if args.samFlagExclude and read.flag & args.samFlagExclude != 0:
            nFiltered += 1
            continue

        tLen = getTLen(read)
        if args.minFragmentLength > 0 and tLen < args.minFragmentLength:
            nFiltered += 1
            continue
        if args.maxFragmentLength > 0 and tLen > args.maxFragmentLength:
            nFiltered += 1
            continue
        if args.ignoreDuplicates:
            # Assuming more or less concordant reads, use the fragment bounds, otherwise the start positions
            if tLen >= 0:
                s = read.pos
                e = s + tLen
            else:
                s = read.pnext
                e = s - tLen
            if read.reference_id != read.next_reference_id:
                e = read.pnext
            if lpos is not None and lpos == read.reference_start \
                    and (s, e, read.next_reference_id, read.is_reverse) in prev_pos:
                nFiltered += 1
                continue
            if lpos != read.reference_start:
                prev_pos.clear()
            lpos = read.reference_start
            prev_pos.add((s, e, read.next_reference_id, read.is_reverse))

        # filterRNAstrand
        if args.filterRNAstrand:
            if read.is_paired:
                if args.filterRNAstrand == 'forward':
                    if read.flag & 144 == 128 or read.flag & 96 == 64:
                        pass
                    else:
                        nFiltered += 1
                        continue
                elif args.filterRNAstrand == 'reverse':
                    if read.flag & 144 == 144 or read.flag & 96 == 96:
                        pass
                    else:
                        nFiltered += 1
                        continue
            else:
                if args.filterRNAstrand == 'forward':
                    if read.flag & 16 == 16:
                        pass
                    else:
                        nFiltered += 1
                        continue
                elif args.filterRNAstrand == 'reverse':
                    if read.flag & 16 == 0:
                        pass
                    else:
                        nFiltered += 1
                        continue

        # Read survived filtering
        ofh.write(read)

    # The results from the workers will get sorted, so get the TID
    tid = fh.get_tid(chrom)

    ofh.close()
    fh.close()
    return tid, start, total, nFiltered, oname


def main(args=None):
    args = parseArguments().parse_args(args)

    bam = openBam(args.bam)
    total = bam.mapped + bam.unmapped
    chrom_sizes = [(x, y) for x, y in zip(bam.header.references, bam.header.lengths)]

    # Filter, writing the results to a bunch of temporary files
    res = mapReduce([args],
                    filterWorker,
                    chrom_sizes,
                    blackListFileName=args.blackListFileName,
                    numberOfProcessors=args.numberOfProcessors,
                    verbose=args.verbose)

    res = sorted(res)  # The temp files are now in order for concatenation
    nFiltered = sum([x[3] for x in res])
    totalSeen = sum([x[2] for x in res])  # The * contig isn't queried

    mode = 'wb'
    if bam.is_cram:
        mode = 'wc'
    obam = pysam.AlignmentFile(args.outFile, mode, template=bam)
    bam.close()

    for tup in res:
        if tup[3] - tup[2] == 0:
            # No alignments, skip
            os.remove(tup[4])
            continue
        bam = pysam.AlignmentFile(tup[4])
        for b in bam.fetch(until_eof=True):
            obam.write(b)
        bam.close()
        os.remove(tup[4])
    obam.close()

    if args.filterMetrics:
        sampleName = args.bam
        if args.label:
            sampleName = args.label
        if args.smartLabels:
            sampleName = smartLabels([args.bam])[0]

        of = open(args.filterMetrics, "w")
        of.write("#bamFilterReads --filterMetrics\n")
        of.write("#File\tReads Remaining\tTotal Initial Reads\n")
        of.write("{}\t{}\t{}\n".format(sampleName, totalSeen - nFiltered, total))
        of.close()

    return 0
