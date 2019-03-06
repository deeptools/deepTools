#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import deeptools.mapReduce as mapReduce
from deeptools import bamHandler
from deeptools import utilities
import sys

debug = 0


def getFractionKept_wrapper(args):
    return getFractionKept_worker(*args)


def getFractionKept_worker(chrom, start, end, bamFile, args, offset):
    """
    Queries the BAM file and counts the number of alignments kept/found in the
    first 50000 bases.
    """
    bam = bamHandler.openBam(bamFile)
    start += offset * 50000
    end = min(end, start + 50000)
    tot = 0
    filtered = 0

    if end <= start:
        return (filtered, tot)

    prev_pos = set()
    lpos = None
    if chrom in bam.references:
        for read in bam.fetch(chrom, start, end):
            tot += 1
            if read.is_unmapped:
                continue

            if args.minMappingQuality and read.mapq < args.minMappingQuality:
                filtered += 1
                continue

            # filter reads based on SAM flag
            if args.samFlagInclude and read.flag & args.samFlagInclude != args.samFlagInclude:
                filtered += 1
                continue
            if args.samFlagExclude and read.flag & args.samFlagExclude != 0:
                filtered += 1
                continue

            # fragment length filtering
            tLen = utilities.getTLen(read)
            if args.minFragmentLength > 0 and tLen < args.minFragmentLength:
                filtered += 1
                continue
            if args.maxFragmentLength > 0 and tLen > args.maxFragmentLength:
                filtered += 1
                continue

            # get rid of duplicate reads that have same position on each of the
            # pairs
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
                    filtered += 1
                    continue
                if lpos != read.reference_start:
                    prev_pos.clear()
                lpos = read.reference_start
                prev_pos.add((s, e, read.next_reference_id, read.is_reverse))

            # If filterRNAstrand is in args, then filter accordingly
            # This is very similar to what's used in the get_fragment_from_read function in the filterRnaStrand class
            if hasattr(args, "filterRNAstrand"):
                if read.is_paired:
                    if args.filterRNAstrand == 'forward':
                        if not ((read.flag & 128 == 128 and read.flag & 16 == 0) or (read.flag & 64 == 64 and read.flag & 32 == 0)):
                            filtered += 1
                            continue
                    elif args.filterRNAstrand == 'reverse':
                        if not (read.flag & 144 == 144 or read.flag & 96 == 96):
                            filtered += 1
                            continue
                else:
                    if args.filterRNAstrand == 'forward' and read.flag & 16 == 0:
                        filtered += 1
                        continue
                    elif args.filterRNAstrand == 'reverse' and read.flag & 16 == 16:
                        filtered += 1
                        continue

    return (filtered, tot)


def fraction_kept(args, stats):
    """
    Count the following:
    (A) The total number of alignments sampled
    (B) The total number of alignments ignored due to any of the following:
        --samFlagInclude
        --samFlagExclude
        --minMappingQuality
        --ignoreDuplicates
        --minFragmentLength
        --maxFragmentLength

    Black list regions are already accounted for. This works by sampling the
    genome (by default, we'll iterate until we sample 1% or 100,000 alignments,
    whichever is smaller (unless there are fewer than 100,000 alignments, in
    which case sample everything).

    The sampling works by dividing the genome into bins and only looking at the
    first 50000 bases. If this doesn't yield sufficient alignments then the bin
    size is halved.
    """
    # Do we even need to proceed?
    if (not args.minMappingQuality or args.minMappingQuality == 0) and \
       (not args.samFlagInclude or args.samFlagInclude == 0) and \
       (not args.samFlagExclude or args.samFlagExclude == 0) and \
       (not args.minFragmentLength or args.minFragmentLength == 0) and \
       (not args.maxFragmentLength or args.maxFragmentLength == 0):
        if hasattr(args, "filterRNAstrand"):
            if args.filterRNAstrand not in ["forward", "reverse"]:
                return 1.0
        else:
            return 1.0

    filtered = 0
    total = 0
    distanceBetweenBins = 2000000
    bam_handle = bamHandler.openBam(args.bam)
    bam_mapped = utilities.bam_total_reads(bam_handle, args.ignoreForNormalization, stats)
    if bam_mapped < 1000000:
        num_needed_to_sample = bam_mapped
    else:
        if 0.1 * bam_mapped >= 1000000:
            num_needed_to_sample = 0.1 * bam_mapped
        else:
            num_needed_to_sample = 1000000
    if args.exactScaling:
        num_needed_to_sample = bam_mapped
    if num_needed_to_sample == bam_mapped:
        distanceBetweenBins = 55000
    if args.ignoreForNormalization:
        chrom_sizes = [(chrom_name, bam_handle.lengths[idx]) for idx, chrom_name in enumerate(bam_handle.references)
                       if chrom_name not in args.ignoreForNormalization]
    else:
        chrom_sizes = list(zip(bam_handle.references, bam_handle.lengths))

    offset = 0
    # Iterate over bins at various non-overlapping offsets until we have enough data
    while total < num_needed_to_sample and offset < np.ceil(distanceBetweenBins / 50000):
        res = mapReduce.mapReduce((bam_handle.filename, args, offset),
                                  getFractionKept_wrapper,
                                  chrom_sizes,
                                  genomeChunkLength=distanceBetweenBins,
                                  blackListFileName=args.blackListFileName,
                                  numberOfProcessors=args.numberOfProcessors,
                                  verbose=args.verbose)

        if len(res):
            foo, bar = np.sum(res, axis=0)
            filtered += foo
            total += bar
        offset += 1

    if total == 0:
        # This should never happen
        total = 1

    return 1.0 - float(filtered) / float(total)


def get_num_kept_reads(args, stats):
    """
    Substracts from the total number of mapped reads in a bamfile
    the proportion of reads that fall into blacklisted regions
    or that are filtered

    :return: integer
    """
    if stats is None:
        bam_handle, mapped, unmapped, stats = bamHandler.openBam(args.bam, returnStats=True, nThreads=args.numberOfProcessors)
    else:
        bam_handle = bamHandler.openBam(args.bam)
    bam_mapped_total = utilities.bam_total_reads(bam_handle, args.ignoreForNormalization, stats)
    if args.blackListFileName:
        blacklisted = utilities.bam_blacklisted_reads(bam_handle, args.ignoreForNormalization,
                                                      args.blackListFileName, args.numberOfProcessors)
        print("There are {0} alignments, of which {1} are completely "
              "within a blacklist region.".format(bam_mapped_total, blacklisted))
        num_kept_reads = bam_mapped_total - blacklisted
    else:
        num_kept_reads = bam_mapped_total
    ftk = fraction_kept(args, stats)
    if ftk < 1:
        num_kept_reads *= ftk
        print("Due to filtering, {0}% of the aforementioned alignments "
              "will be used {1}".format(100 * ftk, num_kept_reads))

    return num_kept_reads, bam_mapped_total


def get_scale_factor(args, stats):
    scale_factor = args.scaleFactor
    bam_mapped, bam_mapped_total = get_num_kept_reads(args, stats)
    if args.normalizeUsing == 'RPGC':
        # Print output, since normalzation stuff isn't printed to stderr otherwise
        sys.stderr.write("normalization: 1x (effective genome size {})\n".format(args.effectiveGenomeSize))

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

        current_coverage = \
            float(bam_mapped * fragment_length) / args.effectiveGenomeSize
        # the scaling sets the coverage to match 1x
        scale_factor *= 1.0 / current_coverage
        if debug:
            print("Estimated current coverage {}".format(current_coverage))
            print("Scaling factor {}".format(args.scaleFactor))

    elif args.normalizeUsing == 'RPKM':
        # Print output, since normalzation stuff isn't printed to stderr otherwise
        sys.stderr.write("normalization: RPKM\n")

        # the RPKM is the # reads per tile / \
        #    ( total reads (in millions) * tile length in Kb)
        million_reads_mapped = float(bam_mapped) / 1e6
        tile_len_in_kb = float(args.binSize) / 1000

        scale_factor *= 1.0 / (million_reads_mapped * tile_len_in_kb)

        if debug:
            print("scale factor using RPKM is {0}".format(args.scaleFactor))

    elif args.normalizeUsing == 'CPM':
        # Print output, since normalzation stuff isn't printed to stderr otherwise
        sys.stderr.write("normalization: CPM\n")

        # the CPM (norm is based on post-filtering total counts of reads in BAM "bam_mapped")
        million_reads_mapped = float(bam_mapped) / 1e6
        scale_factor *= 1.0 / (million_reads_mapped)

        if debug:
            print("scale factor using CPM is {0}".format(args.scaleFactor))

    elif args.normalizeUsing == 'BPM':
        # Print output, since normalzation stuff isn't printed to stderr otherwise
        sys.stderr.write("normalization: BPM\n")
        # the BPM (norm is based on post-filtering total counts of reads in BAM "bam_mapped")
        # sampled_bins_sum = getSampledSum(args.bam)
        tile_len_in_kb = float(args.binSize) / 1000
        tpm_scaleFactor = (bam_mapped / tile_len_in_kb) / 1e6

        scale_factor *= 1 / (tpm_scaleFactor * tile_len_in_kb)
        if debug:
            print("scale factor using BPM is {0}".format(args.scaleFactor))

    else:
        # Print output, since normalzation stuff isn't printed to stderr otherwise
        sys.stderr.write("normalization: none (signal scaled by the fraction of alignments kept after filtering)\n")

        scale_factor *= bam_mapped / float(bam_mapped_total)

    if args.verbose:
        print("Final scaling factor: {}".format(scale_factor))

    return scale_factor
