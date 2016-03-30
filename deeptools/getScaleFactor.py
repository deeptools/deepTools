#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import deeptools.mapReduce as mapReduce
from deeptools import bamHandler
from deeptools import parserCommon

debug = 0


def getFractionKept_wrapper(args):
    return getFractionKept_worker(*args)


def getFractionKept_worker(chrom, start, end, bamFile, args):
    """
    Queries the BAM file and counts the number of alignments kept/found in the
    first 50000 bases.
    """
    bam = bamHandler.openBam(bamFile)
    end = min(end, start + 50000)
    tot = 0
    filtered = 0
    prev_start_pos = None  # to store the start positions
    if chrom in bam.references:
        for read in bam.fetch(chrom, start, end):
            tot += 1
            if args.minMappingQuality and read.mapq < args.minMappingQuality:
                filtered += 1
                continue

            # filter reads based on SAM flag
            if args.samFlagInclude and read.flag & args.samFlagInclude == 0:
                filtered += 1
                continue
            if args.samFlagExclude and read.flag & args.samFlagExclude != 0:
                filtered += 1
                continue

            # get rid of duplicate reads that have same position on each of the
            # pairs
            if args.ignoreDuplicates and prev_start_pos \
                    and prev_start_pos == (read.reference_start, read.pnext, read.is_reverse):
                filtered += 1
                continue
            prev_start_pos = (read.reference_start, read.pnext, read.is_reverse)

    return (filtered, tot)


def fraction_kept(args):
    """
    Count the following:
    (A) The total number of alignments sampled
    (B) The total number of alignments ignored due to any of the following:
        --samFlagInclude
        --samFlagExclude
        --minMappingQuality
        --ignoreDuplicates

    Black list regions are already accounted for. This works by sampling the
    genome (by default, we'll iterate until we sample 1% or 100,000 alignments,
    whichever is smaller (unless there are fewer than 100,000 alignments, in
    which case sample everything).

    The sampling works by dividing the genome into bins and only looking at the
    first 50000 bases. If this doesn't yield sufficient alignments then the bin
    size is halved.
    """
    filtered = 0
    total = 0
    distanceBetweenBins = 2000000
    bam_handle = bamHandler.openBam(args.bam)
    bam_mapped = parserCommon.bam_total_reads(bam_handle, args.ignoreForNormalization)
    num_needed_to_sample = max(bam_mapped if bam_mapped <= 100000 else 0, min(100000, 0.01 * bam_mapped))
    chrom_sizes = zip(bam_handle.references, bam_handle.lengths)

    while total < num_needed_to_sample and distanceBetweenBins > 50000:
        # If we've iterated, then halve distanceBetweenBins
        distanceBetweenBins /= 2
        if distanceBetweenBins < 50000:
            distanceBetweenBins = 50000

        res = mapReduce.mapReduce((bam_handle.filename, args),
                                  getFractionKept_wrapper,
                                  chrom_sizes,
                                  genomeChunkLength=distanceBetweenBins,
                                  blackListFileName=args.blackListFileName,
                                  numberOfProcessors=args.numberOfProcessors,
                                  verbose=args.verbose)

        if len(res):
            filtered, total = np.sum(res, axis=0)

    if total == 0:
        # This should never happen
        total = 1

    return 1.0 - float(filtered) / float(total)


def get_scale_factor(args):
    scale_factor = args.scaleFactor
    bam_handle = bamHandler.openBam(args.bam)
    bam_mapped = parserCommon.bam_total_reads(bam_handle, args.ignoreForNormalization)
    blacklisted = parserCommon.bam_blacklisted_reads(bam_handle, args.ignoreForNormalization, args.blackListFileName)
    if args.verbose:
        print("There are {} alignments, of which {} are completely within a blacklist region.".format(bam_mapped, blacklisted))
    bam_mapped -= blacklisted
    ftk = fraction_kept(args)
    bam_mapped *= ftk
    if args.verbose:
        print("Due to filtering, {}%% of the aforementioned alignments will be used {}".format(100 * ftk, bam_mapped))

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

    if args.verbose:
        print("Final scaling factor: {}".format(scale_factor))

    return scale_factor
