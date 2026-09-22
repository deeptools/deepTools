import os
import shutil
import sys
import tempfile

import numpy as np

# NGS packages
import pyBigWig

# own module
from deeptools import bamHandler, mapReduce
from deeptools.countReadsPerBin import CountReadsPerBin
from deeptools.utilities import getCommonChrNames, toBytes
from deeptools.writeBedGraph import *

old_settings = np.seterr(all='ignore')


def getCoverageFromBigwig(bigwigHandle, chrom, start, end, tileSize,
                          missingDataAsZero=False):
    try:
        coverage = np.asarray(bigwigHandle.values(chrom, start, end))
    except RuntimeError:
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
        skipZeroOverZero=False, missingDataAsZero=False, fixedStep=False):
    r"""
    Writes a bedgraph having as base a number of bam files.

    The given func is called to compute the desired bedgraph value
    using the funcArgs

    tileSize
    """
    if start > end:
        raise NameError(f"start position ({start}) bigger than "
                        f"end position ({end})")

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

    with tempfile.NamedTemporaryFile(delete=False) as _file:
        previousValue = None
        lengthCoverage = len(coverage[0])
        for tileIndex in range(lengthCoverage):

            tileCoverage = []
            for index in range(len(bamOrBwFileList)):
                if smoothLength > 0:
                    vectorStart, vectorEnd = CountReadsPerBin.getSmoothRange(
                        tileIndex, tileSize, smoothLength, lengthCoverage)
                    tileCoverage.append(
                        np.mean(coverage[index][vectorStart:vectorEnd]))
                else:
                    try:
                        tileCoverage.append(coverage[index][tileIndex])
                    except IndexError:
                        sys.exit(f"Chromosome {chrom} probably not in one of the bigwig "
                                 "files. Remove this chromosome from the bigwig file "
                                 "to continue")

            if skipZeroOverZero and np.sum(tileCoverage) == 0:
                previousValue = None
                continue

            value = func(tileCoverage, funcArgs)

            if fixedStep:
                writeStart = start + tileIndex * tileSize
                writeEnd = min(writeStart + tileSize, end)
                try:
                    _file.write(toBytes(f"{chrom}\t{writeStart}\t{writeEnd}\t{value:g}\n"))
                except TypeError:
                    _file.write(toBytes(f"{chrom}\t{writeStart}\t{writeEnd}\t{value}\n"))
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
                            toBytes(f"{chrom}\t{writeStart}\t{writeEnd}\t{previousValue:g}\n"))
                    previousValue = value
                    writeStart = writeEnd
                    writeEnd = min(writeStart + tileSize, end)

        if (not fixedStep) and previousValue and (writeStart != end) and (not np.isnan(previousValue)):
            _file.write(toBytes(f"{chrom}\t{writeStart}\t{end}\t{previousValue:g}\n"))

        tempFileName = _file.name
    return chrom, start, end, tempFileName


def writeBedGraph(
        bamOrBwFileList, outputFileName, fragmentLength,
        func, funcArgs, tileSize=25, region=None, blackListFileName=None, numberOfProcessors=1,
        format="bedgraph", extendPairedEnds=True, missingDataAsZero=False,
        skipZeroOverZero=False, smoothLength=0, fixedStep=False, verbose=False):
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
            bam, mapped, _unmapped, _stats = bamHandler.openBam(indexedFile, returnStats=True, nThreads=numberOfProcessors)
            bamHandles.append(bam)
            mappedList.append(mapped)

    if len(bamHandles):
        genomeChunkLength = getGenomeChunkLength(bamHandles, tileSize, mappedList)
        # check if both bam files correspond to the same species
        # by comparing the chromosome names:
        chromNamesAndSize, __ = getCommonChrNames(bamHandles, verbose=verbose)
    else:
        genomeChunkLength = int(10e6)
        cCommon_number = {}
        chromNamesAndSize = {}
        for fileName, fileFormat in bamOrBwFileList:
            if fileFormat == 'bigwig':
                fh = pyBigWig.open(fileName)
            else:
                continue

            for chromName, size in list(fh.chroms().items()):
                if chromName in chromNamesAndSize:
                    cCommon_number[chromName] += 1
                    if chromNamesAndSize[chromName] != size:
                        print("\nWARNING\n"
                              f"Chromosome {chromName} length reported in the "
                              f"input files differ.\n{chromNamesAndSize[chromName]} for {bamOrBwFileList[0][0]}\n"
                              f"{size} for {fileName}.\n\nThe smallest "
                              "length will be used")
                        chromNamesAndSize[chromName] = min(
                            chromNamesAndSize[chromName], size)
                else:
                    chromNamesAndSize[chromName] = size
                    cCommon_number[chromName] = 1
            fh.close()

        # get the list of common chromosome names and sizes
        if len(bamOrBwFileList) == 1:
            chromNamesAndSize = [(k, v) for k, v in chromNamesAndSize.items()]
        else:
            chromNamesAndSize = [(k, v) for k, v in chromNamesAndSize.items()
                                 if k in cCommon_number and
                                 cCommon_number[k] == len(bamOrBwFileList)]

    if region:
        # in case a region is used, append the tilesize
        region += f":{tileSize}"

    res = mapReduce.mapReduce((tileSize, fragmentLength, bamOrBwFileList,
                               func, funcArgs, extendPairedEnds, smoothLength,
                               skipZeroOverZero, missingDataAsZero, fixedStep),
                              writeBedGraph_wrapper,
                              chromNamesAndSize,
                              genomeChunkLength=genomeChunkLength,
                              region=region,
                              blackListFileName=blackListFileName,
                              numberOfProcessors=numberOfProcessors,
                              verbose=verbose)

    # Determine the sorted order of the temp files
    chrom_order = {}
    for i, _ in enumerate(chromNamesAndSize):
        chrom_order[_[0]] = i
    res = [[chrom_order[x[0]], x[1], x[2], x[3]] for x in res]
    res.sort()

    if format == 'bedgraph':
        with open(outputFileName, 'wb') as of:
            for r in res:
                if r is not None:
                    with open(r[3], 'rb') as _:
                        shutil.copyfileobj(_, of)
                    os.remove(r[3])
    else:
        bedGraphToBigWig(chromNamesAndSize, [x[3] for x in res], outputFileName)
