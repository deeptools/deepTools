#!/usr/bin/env python
import argparse
import pysam
import os
import sys

from deeptools import parserCommon
from deeptools.bamHandler import openBam
from deeptools.mapReduce import mapReduce
from deeptools.utilities import getTLen, smartLabels, getTempFileName
from importlib.metadata import version
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
                           help="A BED or GTF file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant.",
                           metavar="BED file",
                           nargs="+",
                           default="None",
                           required=False)

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


def shiftRead(b, chromDict, args):
    if not b.is_proper_pair:
        return None
    tLen = getTLen(b, notAbs=True)
    start = b.pos
    end = start + b.query_alignment_end
    if b.is_reverse and not b.is_read2:
        end -= args.shift[2]
        deltaTLen = args.shift[3] - args.shift[2]
    elif b.is_reverse and b.is_read2:
        end += args.shift[1]
        deltaTLen = args.shift[1] - args.shift[0]
    elif not b.is_reverse and not b.is_read2:
        start += args.shift[0]
        deltaTLen = args.shift[1] - args.shift[0]
    else:
        start -= args.shift[3]
        deltaTLen = args.shift[3] - args.shift[2]

    # Sanity check
    if end - start < 1:
        if b.is_reverse:
            start = end - 1
        else:
            end = start + 1
    if start < 0:
        start = 0
    if end > chromDict[b.reference_name]:
        end = chromDict[b.reference_name]
    if end - start < 1:
        return None

    # create a new read
    b2 = pysam.AlignedSegment()
    b2.query_name = b.query_name
    b2.flag = b.flag
    b2.reference_id = b.reference_id
    b2.reference_start = start
    b2.mapping_quality = b.mapping_quality
    b2.cigar = ((0, end - start),)  # Returned cigar is only matches
    if tLen < 0:
        b2.template_length = tLen - deltaTLen
    else:
        b2.template_length = tLen + deltaTLen
    b2.next_reference_id = b.next_reference_id
    b2.next_reference_start = b.next_reference_start
    if b.is_proper_pair:
        if b2.is_read2 and b2.is_reverse:
            b2.next_reference_start += args.shift[0]
        elif not b2.is_read2 and b2.is_reverse:
            b2.next_reference_start -= args.shift[3]

    return b2


def filterWorker(arglist):
    chrom, start, end, args, chromDict = arglist
    fh = openBam(args.bam)
    mode = 'wb'
    oname = getTempFileName(suffix='.bam')
    if args.filteredOutReads:
        onameFiltered = getTempFileName(suffix='.bam')
    else:
        onameFiltered = None
    ofh = pysam.AlignmentFile(oname, mode=mode, template=fh)
    if onameFiltered:
        ofiltered = pysam.AlignmentFile(onameFiltered, mode=mode, template=fh)
    else:
        ofiltered = None

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
            if ofiltered:
                ofiltered.write(read)
            continue

        if args.minMappingQuality and read.mapq < args.minMappingQuality:
            nFiltered += 1
            if ofiltered:
                ofiltered.write(read)
            continue

        if args.samFlagInclude and read.flag & args.samFlagInclude != args.samFlagInclude:
            nFiltered += 1
            if ofiltered:
                ofiltered.write(read)
            continue
        if args.samFlagExclude and read.flag & args.samFlagExclude != 0:
            nFiltered += 1
            if ofiltered:
                ofiltered.write(read)
            continue

        tLen = getTLen(read)
        if args.minFragmentLength > 0 and tLen < args.minFragmentLength:
            nFiltered += 1
            if ofiltered:
                ofiltered.write(read)
            continue
        if args.maxFragmentLength > 0 and tLen > args.maxFragmentLength:
            nFiltered += 1
            if ofiltered:
                ofiltered.write(read)
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
                if ofiltered:
                    ofiltered.write(read)
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
                        if ofiltered:
                            ofiltered.write(read)
                        continue
                elif args.filterRNAstrand == 'reverse':
                    if read.flag & 144 == 144 or read.flag & 96 == 96:
                        pass
                    else:
                        nFiltered += 1
                        if ofiltered:
                            ofiltered.write(read)
                        continue
            else:
                if args.filterRNAstrand == 'forward':
                    if read.flag & 16 == 16:
                        pass
                    else:
                        nFiltered += 1
                        if ofiltered:
                            ofiltered.write(read)
                        continue
                elif args.filterRNAstrand == 'reverse':
                    if read.flag & 16 == 0:
                        pass
                    else:
                        nFiltered += 1
                        if ofiltered:
                            ofiltered.write(read)
                        continue

        if args.shift:
            read = shiftRead(read, chromDict, args)
            if not read:
                continue

        # Read survived filtering
        ofh.write(read)

    # The results from the workers will get sorted, so get the TID
    tid = fh.get_tid(chrom)

    ofh.close()
    if ofiltered:
        ofiltered.close()
    fh.close()
    return tid, start, total, nFiltered, oname, onameFiltered


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

    # Remove args:
    # label, smartLabels, genomeChunkLength, ignoreDuplicates.

    print(args)
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
    )
