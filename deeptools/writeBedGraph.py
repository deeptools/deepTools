import os
import sys
import shutil
import numpy as np
import pyBigWig

# own modules
from deeptools import mapReduce
from deeptools.utilities import getCommonChrNames
import deeptools.countReadsPerBin as cr
from deeptools import bamHandler
from deeptools import utilities

debug = 0
old_settings = np.seterr(all='ignore')


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

    Examples
    --------

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
    >>> test_path = os.path.dirname(os.path.abspath(__file__)) + "/test/test_data/"

    >>> outFile = tempfile.NamedTemporaryFile()
    >>> bam_file = test_path +  "testA.bam"

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
    >>> c = WriteBedGraph([bam_file], binLength=bin_length, region=region, stepSize=step_size)
    >>> c.run(function_to_call, funcArgs, outFile.name)
    >>> f = open(outFile.name, 'r')
    >>> f.readlines()
    ['3R\t0\t100\t0\n', '3R\t100\t200\t1.5\n']
    >>> f.close()
    >>> outFile.close()


    """

    def run(self, func_to_call, func_args, out_file_name, blackListFileName=None, format="bedgraph", smoothLength=0):
        r"""
        Given a list of bamfiles, a function and a function arguments,
        this method writes a bedgraph file (or bigwig) file
        for a partition of the genome into tiles of given size
        and a value for each tile that corresponds to the given function
        and that is related to the coverage underlying the tile.

        Parameters
        ----------
        func_to_call : str
            function name to be called to convert the list of coverages computed
            for each bam file at each position into a single value. An example
            is a function that takes the ratio between the coverage of two
            bam files.
        func_args : dict
            dict of arguments to pass to `func`. E.g. {'scaleFactor':1.0}

        out_file_name : str
            name of the file to save the resulting data.

        smoothLength : int
            Distance in bp for smoothing the coverage per tile.


        """
        self.__dict__["smoothLength"] = smoothLength
        getStats = len(self.mappedList) < len(self.bamFilesList)
        bam_handles = []
        for x in self.bamFilesList:
            if getStats:
                bam, mapped, unmapped, stats = bamHandler.openBam(x, returnStats=True, nThreads=self.numberOfProcessors)
                self.mappedList.append(mapped)
                self.statsList.append(stats)
            else:
                bam = bamHandler.openBam(x)
            bam_handles.append(bam)

        genome_chunk_length = getGenomeChunkLength(bam_handles, self.binLength, self.mappedList)
        # check if both bam files correspond to the same species
        # by comparing the chromosome names:
        chrom_names_and_size, non_common = getCommonChrNames(bam_handles, verbose=False)

        if self.region:
            # in case a region is used, append the tilesize
            self.region += ":{}".format(self.binLength)

        for x in list(self.__dict__.keys()):
            if x in ["mappedList", "statsList"]:
                continue
            sys.stderr.write("{}: {}\n".format(x, self.__getattribute__(x)))

        res = mapReduce.mapReduce([func_to_call, func_args],
                                  writeBedGraph_wrapper,
                                  chrom_names_and_size,
                                  self_=self,
                                  genomeChunkLength=genome_chunk_length,
                                  region=self.region,
                                  blackListFileName=blackListFileName,
                                  numberOfProcessors=self.numberOfProcessors)

        # Determine the sorted order of the temp files
        chrom_order = dict()
        for i, _ in enumerate(chrom_names_and_size):
            chrom_order[_[0]] = i
        res = [[chrom_order[x[0]], x[1], x[2], x[3]] for x in res]
        res.sort()

        if format == 'bedgraph':
            out_file = open(out_file_name, 'wb')
            for r in res:
                if r[3]:
                    _foo = open(r[3], 'rb')
                    shutil.copyfileobj(_foo, out_file)
                    _foo.close()
                    os.remove(r[3])
            out_file.close()
        else:
            bedGraphToBigWig(chrom_names_and_size, [x[3] for x in res], out_file_name)

    def writeBedGraph_worker(self, chrom, start, end,
                             func_to_call, func_args,
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
        smoothLength : int
            Distance in bp for smoothing the coverage per tile.
        bed_regions_list: list
            List of tuples of the form (chrom, start, end)
            corresponding to bed regions to be processed.
            If not bed file was passed to the object constructor
            then this list is empty.

        Returns
        -------
        A list of [chromosome, start, end, temporary file], where the temporary file contains the bedgraph results for the region queried.

        Examples
        --------
        >>> test_path = os.path.dirname(os.path.abspath(__file__)) + "/test/test_data/"
        >>> bamFile1 = test_path +  "testA.bam"
        >>> bin_length = 50
        >>> number_of_samples = 0 # overruled by step_size
        >>> func_to_call = scaleCoverage
        >>> funcArgs = {'scaleFactor': 1.0}

        >>> c = WriteBedGraph([bamFile1], bin_length, number_of_samples, stepSize=50)
        >>> tempFile = c.writeBedGraph_worker( '3R', 0, 200, func_to_call, funcArgs)
        >>> f = open(tempFile[3], 'r')
        >>> f.readlines()
        ['3R\t0\t100\t0\n', '3R\t100\t200\t1\n']
        >>> f.close()
        >>> os.remove(tempFile[3])


        """
        if start > end:
            raise NameError("start position ({0}) bigger "
                            "than end position ({1})".format(start, end))

        coverage, _ = self.count_reads_in_region(chrom, start, end)

        _file = open(utilities.getTempFileName(suffix='.bg'), 'w')
        previous_value = None
        line_string = "{}\t{}\t{}\t{:g}\n"
        for tileIndex in range(coverage.shape[0]):

            if self.smoothLength is not None and self.smoothLength > 0:
                vector_start, vector_end = self.getSmoothRange(tileIndex,
                                                               self.binLength,
                                                               self.smoothLength,
                                                               coverage.shape[0])
                tileCoverage = np.mean(coverage[vector_start:vector_end, :], axis=0)
            else:
                tileCoverage = coverage[tileIndex, :]

            value = func_to_call(tileCoverage, func_args)
            """
            # uncomment these lines if fixed step bedgraph is required
            if not np.isnan(value):
                writeStart = start + tileIndex * self.binLength
                writeEnd  =  min(writeStart + self.binLength, end)
                _file.write(line_string.format(chrom, writeStart,
                                               writeEnd, value))
            continue
            """

            if previous_value is None:
                writeStart = start + tileIndex * self.binLength
                writeEnd = min(writeStart + self.binLength, end)
                previous_value = value

            elif previous_value == value:
                writeEnd = min(writeEnd + self.binLength, end)

            elif previous_value != value:
                if not np.isnan(previous_value):
                    _file.write(
                        line_string.format(chrom, writeStart, writeEnd, previous_value))
                previous_value = value
                writeStart = writeEnd
                writeEnd = min(writeStart + self.binLength, end)

        # write remaining value if not a nan
        if previous_value is not None and writeStart != end and not np.isnan(previous_value):
            _file.write(line_string.format(chrom, writeStart,
                                           end, previous_value))

        tempfilename = _file.name
        _file.close()
        return chrom, start, end, tempfilename


def bedGraphToBigWig(chromSizes, bedGraphFiles, bigWigPath):
    """
    Takes a sorted list of bedgraph files and write them to a single bigWig file using pyBigWig.
    The order of bedGraphFiles must match that of chromSizes!
    """
    bw = pyBigWig.open(bigWigPath, "w")
    assert(bw is not None)
    bw.addHeader(chromSizes, maxZooms=10)
    lastChrom = None
    starts = []
    ends = []
    vals = []
    for bg in bedGraphFiles:
        if bg is not None:
            f = open(bg)
            for line in f:
                interval = line.split()
                # Buffer up to a million entries
                if interval[0] != lastChrom or len(starts) == 1000000:
                    if lastChrom is not None:
                        bw.addEntries([lastChrom] * len(starts), starts, ends=ends, values=vals)
                    lastChrom = interval[0]
                    starts = [int(interval[1])]
                    ends = [int(interval[2])]
                    vals = [float(interval[3])]
                else:
                    starts.append(int(interval[1]))
                    ends.append(int(interval[2]))
                    vals.append(float(interval[3]))
            f.close()
            os.remove(bg)
    if len(starts) > 0:
        bw.addEntries([lastChrom] * len(starts), starts, ends=ends, values=vals)
    bw.close()


def getGenomeChunkLength(bamHandles, tile_size, mappedList):
    """
    Tries to estimate the length of the genome sent to the workers
    based on the density of reads per bam file and the number
    of bam files.

    The chunk length should be a multiple of the tileSize

    """

    genomeLength = sum(bamHandles[0].lengths)

    max_reads_per_bp = max([float(x) / genomeLength for x in mappedList])

    # 2e6 is an empirical estimate
    genomeChunkLength = int(min(5e6, int(2e6 / (max_reads_per_bp * len(bamHandles)))))

    genomeChunkLength -= genomeChunkLength % tile_size
    return genomeChunkLength


def scaleCoverage(tile_coverage, args):
    """
    tileCoverage should be an list with only one element
    """
    return args['scaleFactor'] * tile_coverage[0]


def ratio(tile_coverage, args):
    """
    tileCoverage should be an list of two elements
    """
    return float(tile_coverage[0]) / tile_coverage[1]
