import os
import shutil
import numpy as np

# own modules
from deeptools import mapReduce
from deeptools.utilities import getCommonChrNames
from deeptools.countReadsPerBin import getCoverageOfRegion, getSmoothRange
from deeptools import bamHandler
from deeptools import utilities
import config as cfg

debug = 0


def writeBedGraph_wrapper(args):
    return writeBedGraph_worker(*args)


def writeBedGraph_worker(chrom, start, end, tileSize, defaultFragmentLength,
                         bamFilesList, func, funcArgs, extendPairedEnds=True,
                         smoothLength=0, zerosToNans=True,
                         minMappingQuality=None,
                         ignoreDuplicates=False,
                         fragmentFromRead_func=None,
                         centerRead=False):

    r"""
    Writes a bedgraph having as base a number of bam files.

    The given func is called to compute the desired bedgraph value
    using the funcArgs

    tileSize
    >>> test = Tester()
    >>> funcArgs = {'scaleFactor': 1.0}
    >>> tempFile = writeBedGraph_worker( '3R', 0, 200, 50, 0,
    ... [test.bamFile1], scaleCoverage, funcArgs, True, 0, False)
    >>> open(tempFile, 'r').readlines()
    ['3R\t0\t100\t0.00\n', '3R\t100\t200\t1.0\n']
    >>> os.remove(tempFile)

    Test the file being writen for single end reads with
    no extension and no smoothing
    >>> tempFile = writeBedGraph_worker( '3R', 0, 200, 50, 0,
    ... [test.bamFile1], scaleCoverage, funcArgs)
    >>> open(tempFile, 'r').readlines()
    ['3R\t100\t200\t1.0\n']
    >>> os.remove(tempFile)

    Test scaling
    >>> funcArgs = {'scaleFactor': 3.0}
    >>> tempFile = writeBedGraph_worker( '3R', 0, 200, 50, 0,
    ... [test.bamFile1], scaleCoverage, funcArgs)
    >>> open(tempFile, 'r').readlines()
    ['3R\t100\t200\t3.0\n']
    >>> os.remove(tempFile)

    Test ignore duplicates
    >>> funcArgs = {'scaleFactor': 1.0}
    >>> tempFile = writeBedGraph_worker( '3R', 0, 200, 50, 0,
    ... [test.bamFile2], scaleCoverage, funcArgs, ignoreDuplicates=True)
    >>> open(tempFile, 'r').readlines()
    ['3R\t50\t200\t1.0\n']
    >>> os.remove(tempFile)

    Test smoothing
    >>> funcArgs = {'scaleFactor': 1.0}
    >>> tempFile = writeBedGraph_worker( '3R', 100, 200, 20, 0,
    ... [test.bamFile2], scaleCoverage, funcArgs, smoothLength=60)
    >>> open(tempFile, 'r').readlines()
    ['3R\t100\t120\t1.00\n', '3R\t120\t140\t1.67\n', '3R\t140\t160\t2.00\n', '3R\t160\t180\t2.33\n', '3R\t180\t200\t2.0\n']
    >>> os.remove(tempFile)

    Test ratio (needs two bam files)
    >>> funcArgs = {}
    >>> tempFile = writeBedGraph_worker( '3R', 100, 200, 50, 0,
    ... [test.bamFile1, test.bamFile2], ratio , funcArgs)
    >>> open(tempFile, 'r').readlines()
    ['3R\t100\t150\t1.00\n', '3R\t150\t200\t0.5\n']
    >>> os.remove(tempFile)


    Test minMapping quality
    >>> funcArgs = {'scaleFactor': 1.0}
    >>> tempFile = writeBedGraph_worker( '3R', 0, 200, 50, 0,
    ... [test.bamFile2], scaleCoverage, funcArgs, minMappingQuality=40)
    >>> open(tempFile, 'r').readlines()
    ['3R\t150\t200\t1.0\n']
    >>> os.remove(tempFile)

    """
    if start > end:
        raise NameError("start position ({0}) bigger "
                        "than end position ({1})".format(start, end))

    coverage = []
    for bamFile in bamFilesList:
        bamHandle = openBam(bamFile)
        coverage.append(
            getCoverageOfRegion(
                bamHandle, chrom, start, end, tileSize,
                defaultFragmentLength, extendPairedEnds, zerosToNans,
                ignoreDuplicates=ignoreDuplicates,
                minMappingQuality=minMappingQuality,
                fragmentFromRead_func=fragmentFromRead_func,
                centerRead=centerRead))
        bamHandle.close()

    _file = open(utilities.getTempFileName(suffix='.bg'), 'w')
    previousValue = None

    lengthCoverage = len(coverage[0])
    for tileIndex in xrange(lengthCoverage):

        tileCoverage = []
        for index in range(len(bamFilesList)):
            if smoothLength > 0:
                vectorStart, vectorEnd = getSmoothRange(
                    tileIndex, tileSize, smoothLength, lengthCoverage)
                tileCoverage.append(
                    np.mean(coverage[index][vectorStart:vectorEnd]))
            else:
                tileCoverage.append(coverage[index][tileIndex])

        # if zerosToNans == True and sum(tileCoverage) == 0.0:
        #   continue

        value = func(tileCoverage, funcArgs)
        """
        # uncomment this lines if fixed step bedgraph is wanted
        if not  np.isnan(value):
            writeStart = start + tileIndex*tileSize
            writeEnd  =  min(writeStart+tileSize, end)
            _file.write( "%s\t%d\t%d\t%.2f\n" % (chrom, writeStart,
                                                 writeEnd, value) )
        """

        if previousValue is None:
            writeStart = start + tileIndex * tileSize
            writeEnd = min(writeStart + tileSize, end)
            previousValue = value

        elif previousValue == value:
            writeEnd = min(writeEnd + tileSize, end)

        elif previousValue != value:
            if not np.isnan(previousValue):
                _file.write(
                    "{}\t{}\t{}\t{:.2f}\n".format(chrom, writeStart,
                                                  writeEnd, previousValue))
            previousValue = value
            writeStart = writeEnd
            writeEnd = min(writeStart + tileSize, end)

    # write remaining value if not a nan
    if previousValue and writeStart != end and not np.isnan(previousValue):
        _file.write("%s\t%d\t%d\t%.1f\n" % (chrom, writeStart,
                                            end, previousValue))

    tempFileName = _file.name
    _file.close()
    return(tempFileName)


def openBam(bamFile, bamIndex=None):
    return bamHandler.openBam(bamFile, bamIndex)


def bedGraphToBigWig(chromSizes, bedGraphPath, bigWigPath, sort=True):
    """
    takes a bedgraph file, orders it and converts it to
    a bigwig file using command line tools.

    Will fail if the bedGraphToBigWig path changes.
    """

    from tempfile import NamedTemporaryFile
    from os import remove, system

    # destination to save chromosome names and sizes
    _file2 = NamedTemporaryFile(delete=False)

    # bedGraph to bigwig requires the chromosome sizes to be
    # saved into a file
    for chrom, size in chromSizes:
        _file2.write("{}\t{}\n".format(chrom, size))
    _file2.close()

    chrSizesFileName = _file2.name

    # check if the file is empty
    if os.stat(bedGraphPath).st_size < 10:
        import sys
        sys.stderr.write(
            "Error: The generated bedGraphFile was empty. Please adjust\n"
            "your deepTools settings and check your input files.")
        exit(1)

    if sort:
        sort_cmd = cfg.config.get('external_tools', 'sort')
        # temporary file to store sorted bedgraph file
        _file = NamedTemporaryFile(delete=False)
        tempFileName1 = _file.name
        system("{} -k1,1 -k2,2n {} > {}".format(sort_cmd,
                                                bedGraphPath, tempFileName1))
        bedGraphPath = tempFileName1

    bedgraph_to_bigwig = cfg.config.get('external_tools', 'bedgraph_to_bigwig')
    system("{} {} {} {}".format(bedgraph_to_bigwig,
                                bedGraphPath, chrSizesFileName, bigWigPath))

    if sort:
        remove(tempFileName1)

    remove(chrSizesFileName)


def getGenomeChunkLength(bamHandlers, tileSize):
    """
    Tries to estimate the length of the genome sent to the workers
    based on the density of reads per bam file and the number
    of bam files.

    The chunk length should be a multiple of the tileSize

    """

    genomeLength = sum(bamHandlers[0].lengths)

    max_reads_per_bp = max(
        [float(x.mapped) / genomeLength for x in bamHandlers])

    # 2e6 is an empirical estimate
    genomeChunkLength = int(
        min(5e6, int(2e6 / (max_reads_per_bp * len(bamHandlers)))))

    genomeChunkLength -= genomeChunkLength % tileSize
    return genomeChunkLength


def writeBedGraph(bamFilesList, outputFileName, fragmentLength,
                  func, funcArgs, tileSize=25, region=None,
                  numberOfProcessors=None, format="bedgraph",
                  extendPairedEnds=True, zerosToNans=True, smoothLength=0,
                  minMappingQuality=None, ignoreDuplicates=False,
                  fragmentFromRead_func=None,
                  centerRead=False):

    r"""
    Given a list of bamfiles, a function and a function arguments,
    this method writes a bedgraph file (or bigwig) file
    for a partition of the genome into tiles of given size
    and a value for each tile that corresponds to the given function
    and that is related to the coverage underlying the tile.

    >>> test = Tester()
    >>> import tempfile
    >>> outFile = tempfile.NamedTemporaryFile()
    >>> funcArgs = {'scaleFactor': 1.0}
    >>> writeBedGraph( [test.bamFile1], outFile.name,
    ... 0, scaleCoverage, funcArgs, region='3R:0:200')
    >>> open(outFile.name, 'r').readlines()
    ['3R\t100\t200\t1.0\n']
    >>> outFile.close()

    """
    bamHandlers = [openBam(x) for x in bamFilesList]
    genomeChunkLength = getGenomeChunkLength(bamHandlers, tileSize)
    # check if both bam files correspond to the same species
    # by comparing the chromosome names:
    chromNamesAndSize = getCommonChrNames(bamHandlers, verbose=False)

    if region:
        # in case a region is used, append the tilesize
        region += ":{}".format(tileSize)

    res = mapReduce.mapReduce((tileSize, fragmentLength, bamFilesList,
                               func, funcArgs, extendPairedEnds, smoothLength,
                               zerosToNans, minMappingQuality,
                               ignoreDuplicates,
                               fragmentFromRead_func, centerRead),
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


def scaleCoverage(tileCoverage, args):
    """
    tileCoverage should be an list with only one element
    """
    return args['scaleFactor'] * tileCoverage[0]


def ratio(tileCoverage, args):
    """
    tileCoverage should be an list of two elements
    """
    return float(tileCoverage[0]) / tileCoverage[1]


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
        self.root = "./test/test_data/"
        self.bamFile1  = self.root + "testA.bam"
        self.bamFile2  = self.root + "testB.bam"
        self.bamFile_PE  = self.root + "test_paired2.bam"
        self.chrom = '3R'
        global debug
        debug = 0

    def writeBedGraph_worker(self):
        """ prepare arguments for test
        """
        start = 0
        end = 100
        bedGraphStep = 25
        scaleFactors = (1, 1)
        defaultFragmentLength = 10
        return (
            self.chrom, start, end, bedGraphStep,
            scaleFactors, defaultFragmentLength)
