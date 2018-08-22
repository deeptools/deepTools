#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import tempfile
import numpy as np
import sys

# NGS packages
import pyBigWig

# own module
from deeptools import mapReduce
from deeptools.utilities import getCommonChrNames, toBytes
from deeptools.writeBedGraph import *
from deeptools import bamHandler

old_settings = np.seterr(all='ignore')


def getCoverageFromBigwig(bigwigHandle, chrom, start, end, tileSize,
                          missingDataAsZero=False):
    try:
        coverage = np.asarray(bigwigHandle.values(chrom, start, end))
    except TypeError:
        # this error happens when chromosome
        # is not into the bigwig file
        return []
    if coverage is None:
        return []
    if missingDataAsZero is True:
        coverage[np.isnan(coverage)] = 0
    # average the values per bin
    cov = np.array(
        [np.mean(coverage[x:x + tileSize])
         for x in range(0, len(coverage), tileSize)])
    return cov


def writeBedGraph_wrapper(args):
    return writeBedGraph_worker(*args)


def writeBedGraph_worker(
        chrom, start, end, tileSize, defaultFragmentLength,
        bamOrBwFileList, func, funcArgs, extendPairedEnds=True, smoothLength=0,
        missingDataAsZero=False, fixed_step=False):
    r"""
    Writes a bedgraph having as base a number of bam files.

    The given func is called to compute the desired bedgraph value
    using the funcArgs

    tileSize
    """
    if start > end:
        raise NameError("start position ({0}) bigger than "
                        "end position ({1})".format(start, end))

    coverage = []

    for indexFile, fileFormat in bamOrBwFileList:
        if fileFormat == 'bam':
            bamHandle = bamHandler.openBam(indexFile)
            coverage.append(getCoverageFromBam(
                bamHandle, chrom, start, end, tileSize,
                defaultFragmentLength, extendPairedEnds,
                True))
            bamHandle.close()
        elif fileFormat == 'bigwig':
            bigwigHandle = pyBigWig.open(indexFile)
            coverage.append(
                getCoverageFromBigwig(
                    bigwigHandle, chrom, start, end,
                    tileSize, missingDataAsZero))
            bigwigHandle.close()

    # is /dev/shm available?
    # working in this directory speeds the process
    try:
        _file = tempfile.NamedTemporaryFile(dir="/dev/shm", delete=False)
    except OSError:
        _file = tempfile.NamedTemporaryFile(delete=False)

    previousValue = None
    lengthCoverage = len(coverage[0])
    for tileIndex in range(lengthCoverage):

        tileCoverage = []
        for index in range(len(bamOrBwFileList)):
            if smoothLength > 0:
                vectorStart, vectorEnd = getSmoothRange(
                    tileIndex, tileSize, smoothLength, lengthCoverage)
                tileCoverage.append(
                    np.mean(coverage[index][vectorStart:vectorEnd]))
            else:
                try:
                    tileCoverage.append(coverage[index][tileIndex])
                except IndexError:
                    sys.exit("Chromosome {} probably not in one of the bigwig "
                             "files. Remove this chromosome from the bigwig file "
                             "to continue".format(chrom))

#        if  zerosToNans == True and sum(tileCoverage) == 0.0:
#            continue

        value = func(tileCoverage, funcArgs)

        if fixed_step:
            writeStart = start + tileIndex * tileSize
            writeEnd = min(writeStart + tileSize, end)
            try:
                _file.write(toBytes("%s\t%d\t%d\t%.2f\n" % (chrom, writeStart,
                                                            writeEnd, value)))
            except TypeError:
                _file.write(toBytes("{}\t{}\t{}\t{}\n".format(chrom, writeStart,
                                                              writeEnd, value)))
        else:
            if previousValue is None:
                writeStart = start + tileIndex * tileSize
                writeEnd = min(writeStart + tileSize, end)
                previousValue = value

            elif previousValue == value:
                writeEnd = min(writeEnd + tileSize, end)

            elif previousValue != value:
                if not np.isnan(previousValue):
                    _file.write(
                        toBytes("{0}\t{1}\t{2}\t{3:g}\n".format(chrom, writeStart,
                                                                writeEnd, previousValue)))
                previousValue = value
                writeStart = writeEnd
                writeEnd = min(writeStart + tileSize, end)

    if not fixed_step:
        # write remaining value if not a nan
        if previousValue and writeStart != end and \
                not np.isnan(previousValue):
            _file.write(toBytes("{0}\t{1}\t{2}\t{3:g}\n".format(chrom, writeStart,
                                                                end, previousValue)))

    tempFileName = _file.name
    _file.close()
    return chrom, start, end, tempFileName


def writeBedGraph(
        bamOrBwFileList, outputFileName, fragmentLength,
        func, funcArgs, tileSize=25, region=None, blackListFileName=None, numberOfProcessors=1,
        format="bedgraph", extendPairedEnds=True, missingDataAsZero=False,
        smoothLength=0, fixed_step=False, verbose=False):
    r"""
    Given a list of bamfiles, a function and a function arguments,
    this method writes a bedgraph file (or bigwig) file
    for a partition of the genome into tiles of given size
    and a value for each tile that corresponds to the given function
    and that is related to the coverage underlying the tile.

    """
    bamHandles = []
    mappedList = []
    for indexedFile, fileFormat in bamOrBwFileList:
        if fileFormat == 'bam':
            bam, mapped, unmapped, stats = bamHandler.openBam(indexedFile, returnStats=True, nThreads=numberOfProcessors)
            bamHandles.append(bam)
            mappedList.append(mapped)

    if len(bamHandles):
        genomeChunkLength = getGenomeChunkLength(bamHandles, tileSize, mappedList)
        # check if both bam files correspond to the same species
        # by comparing the chromosome names:
        chromNamesAndSize, __ = getCommonChrNames(bamHandles, verbose=verbose)
    else:
        genomeChunkLength = int(10e6)
        cCommon = []
        chromNamesAndSize = {}
        for fileName, fileFormat in bamOrBwFileList:
            if fileFormat == 'bigwig':
                fh = pyBigWig.open(fileName)
            else:
                continue

            for chromName, size in list(fh.chroms().items()):
                if chromName in chromNamesAndSize:
                    cCommon.append(chromName)
                    if chromNamesAndSize[chromName] != size:
                        print("\nWARNING\n"
                              "Chromosome {} length reported in the "
                              "input files differ.\n{} for {}\n"
                              "{} for {}.\n\nThe smallest "
                              "length will be used".format(
                                  chromName, chromNamesAndSize[chromName],
                                  bamOrBwFileList[0][0], size, fileName))
                        chromNamesAndSize[chromName] = min(
                            chromNamesAndSize[chromName], size)
                else:
                    chromNamesAndSize[chromName] = size
            fh.close()

        # get the list of common chromosome names and sizes
        chromNamesAndSize = [(k, v) for k, v in chromNamesAndSize.items()
                             if k in cCommon]

    if region:
        # in case a region is used, append the tilesize
        region += ":{}".format(tileSize)

    res = mapReduce.mapReduce((tileSize, fragmentLength, bamOrBwFileList,
                               func, funcArgs, extendPairedEnds, smoothLength,
                               missingDataAsZero, fixed_step),
                              writeBedGraph_wrapper,
                              chromNamesAndSize,
                              genomeChunkLength=genomeChunkLength,
                              region=region,
                              blackListFileName=blackListFileName,
                              numberOfProcessors=numberOfProcessors,
                              verbose=verbose)

    # Determine the sorted order of the temp files
    chrom_order = dict()
    for i, _ in enumerate(chromNamesAndSize):
        chrom_order[_[0]] = i
    res = [[chrom_order[x[0]], x[1], x[2], x[3]] for x in res]
    res.sort()

    if format == 'bedgraph':
        of = open(outputFileName, 'wb')
        for r in res:
            if r is not None:
                _ = open(r[3], 'rb')
                shutil.copyfileobj(_, of)
                _.close()
                os.remove(r[3])
        of.close()
    else:
        bedGraphToBigWig(chromNamesAndSize, [x[3] for x in res], outputFileName)
