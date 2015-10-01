import os
import shutil
import numpy as np

# own modules
from deeptools import mapReduce
from deeptools.utilities import getCommonChrNames
import deeptools.countReadsPerBin as cr
from deeptools import bamHandler
from deeptools import utilities
import config as cfg

debug = 0

def writeBedGraph_wrapper(args):
    """
    Passes the arguments to writeBedGraph_worker.
    This is a step required given
    the constrains from the multiprocessing module.
    The args var, contains as first element the 'self' value
    from the WriteBedGraph object

    """
    return WriteBedGraph.writeBedGraph_worker(*args)

class WriteBedGraph(cr.CountReadsPerBin):
    r"""Reads bam files coverages and writes a bedgraph or bigwig file

    Extends the CountReadsPerBin object such that the coverage
    of bam files is writen to multiple bedgraph files at once.

    The bedgraph files are later merge into one and converted
    into a bigwig file if necessary.

    The constructor arguments are the same as for CountReadsPerBin. However,
    when calling the `run` method, the following parameters have
    to be passed

    Example
    _______

    Given the following distribution of reads that cover 200 on
    a chromosome named '3R'::

          0                              100                           200
          |------------------------------------------------------------|
        A                                ===============
                                                        ===============


        B                 ===============               ===============
                                         ===============
                                                        ===============

    >>> import tempfile
    >>> outFile = tempfile.NamedTemporaryFile()
    >>> bam_file =  "./test/test_data/testA.bam"

    For the example a simple scaling function is going to be used. This function
    takes the coverage found at each region and multiplies it to the scaling factor.
    In this case the scaling factor is 1.5

    >>> function_to_call = scaleCoverage
    >>> funcArgs = {'scaleFactor': 1.5}

    Restrict process to a region between positions 0 and 200 of chromosome 3R

    >>> region = '3R:0:200'

    Set up such that coverage is computed for consecutive bins of length 25 bp
    >>> bin_length = 25
    >>> step_size = 25

    >>> num_sample_sites = 0 #overruled by step_size
    >>> fragment_length = 0 # if less thatn read length, then read length will be used instead

    >>> c = WriteBedGraph([bam_file], binLength=bin_length, defaultFragmentLength=fragment_length,
    ... region=region, stepSize=step_size)
    >>> c.run(function_to_call, funcArgs, outFile.name)
    >>> open(outFile.name, 'r').readlines()
    ['3R\t0\t100\t0.00\n', '3R\t100\t200\t1.5\n']
    >>> outFile.close()


    """

    def run(self, func_to_call, func_args, out_file_name, format="bedgraph", smooth_length=0):

        r"""
        Given a list of bamfiles, a function and a function arguments,
        this method writes a bedgraph file (or bigwig) file
        for a partition of the genome into tiles of given size
        and a value for each tile that corresponds to the given function
        and that is related to the coverage underlying the tile.

       Parameters
        ----------
        chrom : str
            Chrom name
        start : int
            start coordinate
        end : int
            end coordinate
        func_to_call : str
            function name to be called to convert the list of coverages computed
            for each bam file at each position into a single value. An example
            is a function that takes the ratio between the coverage of two
            bam files.
        func_args : dict
            dict of arguments to pass to `func`. E.g. {'scaleFactor':1.0}

        out_file_name : str
            name of the file to save the resulting data.

        smooth_length : int
            Distance in bp for smoothing the coverage per tile.


        """
        bamHandlers = [bamHandler.openBam(x) for x in self.bamFilesList]
        genomeChunkLength = getGenomeChunkLength(bamHandlers, self.binLength)
        # check if both bam files correspond to the same species
        # by comparing the chromosome names:
        chromNamesAndSize = getCommonChrNames(bamHandlers, verbose=False)

        if self.region:
            # in case a region is used, append the tilesize
            self.region += ":{}".format(self.binLength)

        for x in self.__dict__.keys():
            print "{}: {}".format(x, self.__getattribute__(x))

        res = mapReduce.mapReduce([func_to_call, func_args],
                                  writeBedGraph_wrapper,
                                  chromNamesAndSize,
                                  self_=self,
                                  genomeChunkLength=genomeChunkLength,
                                  region=self.region,
                                  numberOfProcessors=self.numberOfProcessors)

        # concatenate intermediary bedgraph files
        outFile = open(out_file_name + ".bg", 'wb')
        for tempFileName in res:
            if tempFileName:
                # concatenate all intermediate tempfiles into one
                # bedgraph file
                shutil.copyfileobj(open(tempFileName, 'rb'), outFile)
                os.remove(tempFileName)

        bedGraphFile = outFile.name
        outFile.close()
        if format == 'bedgraph':
            os.rename(bedGraphFile, out_file_name)
            if debug:
                print "output file: %s" % (out_file_name)
        else:
            bedGraphToBigWig(
                chromNamesAndSize, bedGraphFile, out_file_name, True)
            if debug:
                print "output file: %s" % (out_file_name)
            os.remove(bedGraphFile)


    def writeBedGraph_worker(self, chrom, start, end,
                             func_to_call, func_args, smooth_length=0,
                             bed_regions_list=None):

        r"""Writes a bedgraph based on the read coverage found on bamFiles

        The given func is called to compute the desired bedgraph value
        using the funcArgs

        Parameters
        ----------
        chrom : str
            Chrom name
        start : int
            start coordinate
        end : int
            end coordinate
        func_to_call : str
            function name to be called to convert the list of coverages computed
            for each bam file at each position into a single value. An example
            is a function that takes the ratio between the coverage of two
            bam files.
        func_args : dict
            dict of arguments to pass to `func`.
        smooth_length : int
            Distance in bp for smoothing the coverage per tile.
        bed_regions_list: list
            List of tuples of the form (chrom, start, end)
            corresponding to bed regions to be processed.
            If not bed file was passed to the object constructor
            then this list is empty.

        Returns
        -------
        temporary file with the bedgraph results for the region queried.

        Example
        -------
        >>> bamFile1  = "./test/test_data/testA.bam"
        >>> bin_length = 50
        >>> number_of_samples = 0 # overruled by step_size
        >>> default_fragment_length = 0 # if < read length, then read length is used instead
        >>> func_to_call = scaleCoverage
        >>> funcArgs = {'scaleFactor': 1.0}

        >>> c = WriteBedGraph([bamFile1], bin_length, number_of_samples,
        ... default_fragment_length, stepSize=50, skipZeros=False)
        >>> tempFile = c.writeBedGraph_worker( '3R', 0, 200, func_to_call, funcArgs)
        >>> open(tempFile, 'r').readlines()
        ['3R\t0\t100\t0.00\n', '3R\t100\t200\t1.0\n']
        >>> os.remove(tempFile)


        """
        if start > end:
            raise NameError("start position ({0}) bigger "
                            "than end position ({1})".format(start, end))

        coverage = []
        bamHandlers = [bamHandler.openBam(bam) for bam in self.bamFilesList]
        for bam in bamHandlers:
            coverage.append(
                self.get_coverage_of_region(bam, chrom, start, end, self.binLength))
            bam.close()

        _file = open(utilities.getTempFileName(suffix='.bg'), 'w')
        previousValue = None

        lengthCoverage = len(coverage[0])
        for tileIndex in xrange(lengthCoverage):

            tileCoverage = []
            for index in range(len(self.bamFilesList)):
                if smooth_length > 0:
                    vectorStart, vectorEnd = self.getSmoothRange(tileIndex,
                                                                self.binLength,
                                                                smooth_length,
                                                                lengthCoverage)
                    tileCoverage.append(
                        np.mean(coverage[index][vectorStart:vectorEnd]))
                else:
                    tileCoverage.append(coverage[index][tileIndex])

            # if zerosToNans == True and sum(tileCoverage) == 0.0:
            #   continue

            value = func_to_call(tileCoverage, func_args)
            """
            # uncomment this lines if fixed step bedgraph is wanted
            if not  np.isnan(value):
                writeStart = start + tileIndex*self.binLength
                writeEnd  =  min(writeStart+self.binLength, end)
                _file.write( "%s\t%d\t%d\t%.2f\n" % (chrom, writeStart,
                                                     writeEnd, value) )
            """

            if previousValue is None:
                writeStart = start + tileIndex * self.binLength
                writeEnd = min(writeStart + self.binLength, end)
                previousValue = value

            elif previousValue == value:
                writeEnd = min(writeEnd + self.binLength, end)

            elif previousValue != value:
                if not np.isnan(previousValue):
                    _file.write(
                        "{}\t{}\t{}\t{:.2f}\n".format(chrom, writeStart,
                                                      writeEnd, previousValue))
                previousValue = value
                writeStart = writeEnd
                writeEnd = min(writeStart + self.binLength, end)

        # write remaining value if not a nan
        if previousValue and writeStart != end and not np.isnan(previousValue):
            _file.write("%s\t%d\t%d\t%.1f\n" % (chrom, writeStart,
                                                end, previousValue))

        tempFileName = _file.name
        _file.close()
        return(tempFileName)


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
        system("LC_ALL=C {} -k1,1 -k2,2n {} > {}".format(sort_cmd,
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
