#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os
import shutil
import tempfile
import numpy as np

# NGS packages
from bx.bbi.bigwig_file import BigWigFile

# own module
import mapReduce
from utilities import getCommonChrNames
from writeBedGraph import *
import config as cfg


def getCoverageFromBam(bamHandle, chrom, start, end, tileSize,
                       defaultFragmentLength, extendPairedEnds=True,
                       zerosToNans=True):
    return getCoverageOfRegion(bamHandle, chrom, start, end, tileSize,
                               defaultFragmentLength,
                               extendPairedEnds=extendPairedEnds,
                               zerosToNans=zerosToNans)


def getCoverageFromBigwig(bigwigHandle, chrom, start, end, tileSize,
                          missingDataAsZero=False):
    try:
        coverage = bigwigHandle.get_as_array(chrom, start, end)
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
    >>> test = Tester()
    >>> funcArgs = {'scaleFactor': 1.0}
    >>> tempFile = writeBedGraph_worker(
    ... '3R', 0, 200, 50, 0, [(test.bamFile1,'bam')],
    ... scaleCoverage, funcArgs, True, 0, False)
    >>> open(tempFile, 'r').readlines()
    ['3R\t100\t200\t1.0\n']
    >>> os.remove(tempFile)

    Test the file being writen for single end reads with no
    extension and no smoothing
    >>> tempFile = writeBedGraph_worker(
    ... '3R', 0, 200, 50, 0, [(test.bamFile1,'bam')],
    ... scaleCoverage, funcArgs)
    >>> open(tempFile, 'r').readlines()
    ['3R\t100\t200\t1.0\n']
    >>> os.remove(tempFile)

    Test scaling
    >>> funcArgs = {'scaleFactor': 3.0}
    >>> tempFile = writeBedGraph_worker(
    ... '3R', 0, 200, 50, 0, [(test.bamFile1,'bam')],
    ... scaleCoverage, funcArgs)
    >>> open(tempFile, 'r').readlines()
    ['3R\t100\t200\t3.0\n']
    >>> os.remove(tempFile)

    Test smoothing
    >>> funcArgs = {'scaleFactor': 1.0}
    >>> tempFile = writeBedGraph_worker(
    ... '3R', 100, 200, 20, 0, [(test.bamFile2,'bam')],
    ... scaleCoverage, funcArgs, smoothLength=60)
    >>> open(tempFile, 'r').readlines()
    ['3R\t100\t120\t1.00\n', '3R\t120\t140\t1.67\n', '3R\t140\t160\t2.00\n', '3R\t160\t180\t2.33\n', '3R\t180\t200\t2.0\n']
    >>> os.remove(tempFile)

    Test ratio (needs two bam files)
    >>> funcArgs = {}
    >>> tempFile = writeBedGraph_worker(
    ... '3R', 100, 200, 50, 0, [(test.bamFile1, 'bam'),
    ... (test.bamFile2, 'bam')], ratio , funcArgs)
    >>> open(tempFile, 'r').readlines()
    ['3R\t100\t150\t1.00\n', '3R\t150\t200\t0.5\n']
    >>> os.remove(tempFile)
    """
    if start > end:
        raise NameError("start position ({0}) bigger than "
                        "end position ({1})".format(start, end))

    coverage = []

    for indexFile, fileFormat in bamOrBwFileList:
        if fileFormat == 'bam':
            bamHandle = openBam(indexFile)
            coverage.append(getCoverageFromBam(
                bamHandle, chrom, start, end, tileSize,
                defaultFragmentLength, extendPairedEnds,
                True))
            bamHandle.close()
        elif fileFormat == 'bigwig':
            fileHandle = open(indexFile, 'r')
            bigwigHandle = BigWigFile(file=fileHandle)
            coverage.append(
                getCoverageFromBigwig(
                    bigwigHandle, chrom, start, end,
                    tileSize, missingDataAsZero))
            fileHandle.close()

    # is /dev/shm available?
    # working in this directory speeds the process
    try:
        _file = tempfile.NamedTemporaryFile(dir="/dev/shm", delete=False)
    except OSError:
        _file = tempfile.NamedTemporaryFile(delete=False)

    previousValue = None
    lengthCoverage = len(coverage[0])
    for tileIndex in xrange(lengthCoverage):

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
                    print "Chromosome {} probably not in one of the bigwig " \
                        "files. Remove this chromosome from the bigwig file " \
                        "to continue".format(chrom)
                    exit(0)

#        if  zerosToNans == True and sum(tileCoverage) == 0.0:
#            continue

        value = func(tileCoverage, funcArgs)

        if fixed_step:
            writeStart = start + tileIndex*tileSize
            writeEnd = min(writeStart + tileSize, end)
            try:
                _file.write("%s\t%d\t%d\t%.2f\n" % (chrom, writeStart,
                                                    writeEnd, value))
            except TypeError:
                _file.write("{}\t{}\t{}\t{}\n".format(chrom, writeStart,
                                                      writeEnd, value))
        else:
            if previousValue is None:
                writeStart = start + tileIndex*tileSize
                writeEnd = min(writeStart + tileSize, end)
                previousValue = value

            elif previousValue == value:
                writeEnd = min(writeEnd + tileSize, end)

            elif previousValue != value:
                if not np.isnan(previousValue):
                    _file.write(
                        "%s\t%d\t%d\t%.2f\n" % (chrom, writeStart,
                                                writeEnd, previousValue))
                previousValue = value
                writeStart = writeEnd
                writeEnd = min(writeStart + tileSize, end)

    if not fixed_step:
        # write remaining value if not a nan
        if previousValue and writeStart != end and \
                not np.isnan(previousValue):
            _file.write("%s\t%d\t%d\t%.1f\n" % (chrom, writeStart,
                                                end, previousValue))

#        """
    tempFileName = _file.name
    _file.close()
    return(tempFileName)


def writeBedGraph(
        bamOrBwFileList, outputFileName, fragmentLength,
        func, funcArgs, tileSize=25, region=None, numberOfProcessors=None,
        format="bedgraph", extendPairedEnds=True, missingDataAsZero=False,
        smoothLength=0, fixed_step=False):

    r"""
    Given a list of bamfiles, a function and a function arguments,
    this method writes a bedgraph file (or bigwig) file
    for a partition of the genome into tiles of given size
    and a value for each tile that corresponds to the given function
    and that is related to the coverage underlying the tile.

    >>> test = Tester()
    >>> outFile = tempfile.NamedTemporaryFile()
    >>> funcArgs = {'scaleFactor': 1.0}
    >>> writeBedGraph([(test.bamFile1, 'bam')], outFile.name,
    ... 0, scaleCoverage, funcArgs, region='3R:0:200')
    >>> open(outFile.name, 'r').readlines()
    ['3R\t100\t200\t1.0\n']
    >>> outFile.close()


    """

    bigwig_info = cfg.config.get('external_tools', 'bigwig_info')
    bamHandlers = [openBam(indexedFile) for
                   indexedFile,
                   fileFormat in bamOrBwFileList if fileFormat == 'bam']
    if len(bamHandlers):
        genomeChunkLength = getGenomeChunkLength(bamHandlers, tileSize)
        # check if both bam files correspond to the same species
        # by comparing the chromosome names:
        chromNamesAndSize = getCommonChrNames(bamHandlers, verbose=False)
    else:
        genomeChunkLength = int(10e6)
        bigwigs = [fileName for fileName,
                   fileFormat in bamOrBwFileList if fileFormat == 'bigwig']
        cCommon = []
        chromNamesAndSize = {}
        for bw in bigwigs:
            inBlock = False
            for line in os.popen(
                    "{} -chroms {}".format(bigwig_info, bw)).readlines():

                if line[0:10] == "chromCount":
                    inBlock = True
                    continue

                if line[0:5] == "bases":
                    break
                if inBlock:
                    chromName, id, size = line.strip().split(" ")
                    size = int(size)
                    if chromName in chromNamesAndSize:
                        cCommon.append(chromName)
                        if chromNamesAndSize[chromName] != size:
                            print "\nWARNING\n" \
                                "Chromosome {} length reported in the " \
                                "bigwig files differ.\n{} for {}\n" \
                                "{} for {}.\n\nThe smallest " \
                                "length will be used".format(
                                chromName, chromNamesAndSize[chromName],
                                bigwigs[0], size, bigwigs[1])
                            chromNamesAndSize[chromName] = min(
                                chromNamesAndSize[chromName], size)
                    else:
                        chromNamesAndSize[chromName] = size

        # get the list of common chromosome names and sizes
        chromNamesAndSize = [(k, v) for k, v in chromNamesAndSize.iteritems()
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
                              numberOfProcessors=numberOfProcessors)

    # concatenate intermediary bedgraph files
    outFile = open(outputFileName + ".bg", 'wb')
    for tempFileName in res:
        if tempFileName:
            # concatenate all intermediate tempfiles into one
            # bedgraph file
            shutil.copyfileobj(open(tempFileName, 'rb'), outFile)
            os.remove(tempFileName)

    bedGraphFile = outFile.name
    outFile.close()
    if format == 'bedgraph':
        os.rename(bedGraphFile, outputFileName)
        if debug:
            print "output file: %s" % (outputFileName)
    else:
        bedGraphToBigWig(
            chromNamesAndSize, bedGraphFile, outputFileName, False)
        if debug:
            print "output file: %s" % (outputFileName)
        os.remove(bedGraphFile)


class Tester():
    def __init__(self):
        """
        The distribution of reads between the two bam files is as follows.

        They cover 200 bp

          0                              100                           200
          |------------------------------------------------------------|
        A                                ===============
                                                        ===============


        B                 ===============               ===============                           
                                         ===============
                                                        ===============
        """
        self.root = "/data/projects/ramirez/tools/deepTools/deeptools/test/test_data/"
        self.bamFile1  = self.root + "testA.bam"
        self.bamFile2  = self.root + "testB.bam"
        self.bamFile_PE  = self.root + "test_paired2.bam"
        self.chrom = '3R'
        global debug
        debug = 0

    def getRead(self, readType):
        """ prepare arguments for test
        """
        bam = openBam(self.bamFile_PE)
        if readType == 'paired-reverse':
            read = [x for x in bam.fetch('chr2', 5000081, 5000082 )][0]
        elif readType == 'single-forward':
            read = [x for x in bam.fetch('chr2', 5001491, 5001492 )][0]
        elif readType == 'single-reverse':
            read = [x for x in bam.fetch('chr2', 5001700, 5001701 )][0]
        else:  # by default a forward paired read is returned
            read = [x for x in bam.fetch('chr2', 5000027, 5000028 )][0]
        return read

    def writeBedGraph_worker(self):
        """ prepare arguments for test
        """
        start = 0
        end = 100
        bedGraphStep = 25
        scaleFactors = (1, 1)
        defaultFragmentLength = 10
        return (self.chrom, start, end, bedGraphStep, scaleFactors, defaultFragmentLength)
