import numpy as np
import time, sys, multiprocessing

# deepTools packages
from utilities import getCommonChrNames
import bamHandler
import mapReduce

debug = 0




def countReadsInRegions_wrapper(args):
    """
    Passes the arguments to countReadsInRegions_worker.
    This is a step required given
    the constrains from the multiprocessing module.
    The args var, contains as first element the 'self' value
    from the countReadsPerBin object

    """
    return CountReadsPerBin.count_reads_in_region(*args)


class CountReadsPerBin(object):
    r"""Collects coverage over multiple bam files using multiprocessing

    This function collects read counts (coverage) from several bam files and returns
    an numpy array with the results. This class uses multiprocessing to compute the coverage.

    Parameters
    ----------
    bamFilesList : list
        List containing the names of indexed bam files. E.g. ['file1.bam', 'file2.bam']

    binLength : int
        Length of the window/bin. This value is overruled by ``bedFile`` if present.

    numberOfSamples : int
        Total number of samples. The genome is divided into ``numberOfSamples``, each
        with a window/bin length equal to ``binLength``. This value is overruled
        by ``stepSize`` in case such value is present and by ``bedFile`` in which
        case the number of samples and bins are defined in the bed file

    defaultFragmentLength : int
        fragment length to extend reads that are not paired. Paired reads are extended to
        the fragment length defined by the mate distance. For Illumina reads, usual values
        are around 300. This value can be determined using the peak caller MACS2 or can be
        approximated by the fragment lengths computed when preparing the library for sequencing.

    numberOfProcessors : int
        Number of processors to use. Default is 4

    verbose : bool
        Output messages. Default: False

    region : str
        Region to limit the computation in the form chrom:start:end.

    bedFile : file_handle
        File handle of a bed file containing the regions for wich to compute the coverage. This option
        overrules ``binLength``, ``numberOfSamples`` and ``stepSize``.
    extendPairedEnds : bool
        Whether coverage should be computed for the extended read length (i.e. the region covered
        by the two mates or the regions expected to be covered by single-reads). Default: true

    minMappingQuality : int
        Reads of a mapping quality less than the give value are not considered. Default: None

    ignoreDuplicates : bool
        Whether read duplicates (same start, end position. If paired-end, same start-end for mates) are
        to be excluded. Default: false

    chrToSkip: list
        List with names of chromosomes that do not want to be included in the coverage computation.
        This is useful to remove unwanted chromosomes (e.g. 'random' or 'Het').

    stepSize : int
        the positions for which the coverage is computed are defined as follows:
        ``range(start, end, stepSize)``. Thus, a stepSize of 1, will compute
        the coverage at each base pair. If the stepSize is equal to the
        binLength then the coverage is computed for consecutive bins. If seepSize is
        smaller than the binLength, then teh bins will overlap.

    center_read : bool
        Determines if reads should be centered with respect to the fragment length.

    samFlag_include : int
        Extracts only those reads having the SAM flag. For example, to get only
        reads that are the first mates a samFlag of 64 could be used. Similarly, the
        samFlag_include can be used to select only reads mapping on the reverse strand
        or to get only properly paired reads.

    samFlag_exclude : int
        Removes reads that match the SAM flag. For example to get all reads
        that map to the forward strand a samFlag_exlude 16 should be used. Which
        translates into exclude all reads that map to the reverse strand.

    zerosToNans : bool
        If true, zero values encountered are transformed to Nans. Default false.


    Returns
    -------
    numpy array

        Each row correspond to each bin/bed region and each column correspond to each of
        the bamFiles.


    Examples
    --------

    The test data contains reads for 200 bp.

    >>> test = Tester()

    The transpose function is used to get a nicer looking output.
    The first line corresponds to the number of reads per bin in bam file 1

    >>> c = CountReadsPerBin([test.bamFile1, test.bamFile2],
    ... 50, 4, 0)
    >>> np.transpose(c.run())
    array([[ 0.,  0.,  1.,  1.],
           [ 0.,  1.,  1.,  2.]])
    """

    def __init__(self, bamFilesList, binLength=50, numberOfSamples=None,
                 defaultFragmentLength=300, numberOfProcessors=1,
                 verbose=False, region=None,
                 bedFile=None, extendPairedEnds=True,
                 minMappingQuality=None,
                 ignoreDuplicates=False,
                 chrsToSkip=[],
                 stepSize=None,
                 center_read=False,
                 samFlag_include=None,
                 samFlag_exclude=None,
                 zerosToNans=False,
                 smoothLength=0):


        self.bamFilesList = bamFilesList
        self.binLength = binLength
        self.numberOfSamples = numberOfSamples
        self.defaultFragmentLength = defaultFragmentLength
        self.numberOfProcessors = numberOfProcessors
        self.verbose = verbose
        self.region = region
        self.bedFile = bedFile
        self.extendPairedEnds = extendPairedEnds
        self.minMappingQuality = minMappingQuality
        self.ignoreDuplicates = ignoreDuplicates
        self.chrsToSkip =  chrsToSkip
        self.stepSize = stepSize
        self.center_read = center_read
        self.samFlag_include = samFlag_include
        self.samFlag_exclude = samFlag_exclude
        self.zerosToNans = zerosToNans
        self.smoothLength=smoothLength

        # check that wither numberOfSamples or stepSize are set
        if numberOfSamples is None and stepSize is None and bedFile is None:
            raise ValueError("either stepSize, numberOfSamples or beFile has to be set")

        self.maxPairedFragmentLength = 2 * self.defaultFragmentLength if \
                                    self.defaultFragmentLength > 100 else 1000

    def run(self):
        # Try to determine an optimal fraction of the genome (chunkSize) that is sent to
        # workers for analysis. If too short, too much time is spend loading the files
        # if too long, some processors end up free.
        # the following values are empirical

        bamFilesHandlers = [bamHandler.openBam(x) for x in self.bamFilesList]
        chromSizes = getCommonChrNames(bamFilesHandlers, verbose=self.verbose)

        # skip chromosome in the list. This is usually for the
        # X chromosome which may have either one copy  in a male sample
        # or a mixture of male/female and is unreliable.
        # Also the skip may contain heterochromatic regions and
        # mitochondrial DNA
        if len(self.chrsToSkip): chromSizes = [x for x in chromSizes if x[0] not in self.chrsToSkip]

        chrNames, chrLengths = zip(*chromSizes)

        genomeSize = sum(chrLengths)
        if self.stepSize is None:
            if self.region is None:
                self.stepSize = max(int(float(genomeSize) / self.numberOfSamples), 1)
            else:
                # compute the step size, based on the number of samples
                # and the length of the region studied
                (chrom, start, end) = mapReduce.getUserRegion(chromSizes, self.region)[:3]
                self.stepSize = max(int(float(end-start) / self.numberOfSamples), 1)

        # number of samples is better if large
        if np.mean(chrLengths) < self.stepSize:
            min_num_of_samples = int(genomeSize / np.mean(chrLengths))
            raise ValueError("numberOfSamples has to be bigger than {} ".format(min_num_of_samples))

        max_mapped = max([x.mapped for x in bamFilesHandlers])

        reads_per_bp = float(max_mapped) / genomeSize
        # chunkSize =  int(100 / ( reads_per_bp  * len(bamFilesList)) )


        chunkSize = int(self.stepSize * 1e3 / (reads_per_bp * len(bamFilesHandlers)))
        [bam_h.close() for bam_h in bamFilesHandlers]

        if self.verbose:
            print "step size is {}".format(self.stepSize)

        if self.region:
            # in case a region is used, append the tilesize
            self.region += ":{}".format(self.binLength)

        # use map reduce to call countReadsInRegions_wrapper
        imap_res = mapReduce.mapReduce([],
                                       countReadsInRegions_wrapper,
                                       chromSizes,
                                       self_=self,
                                       genomeChunkLength=chunkSize,
                                       bedFile=self.bedFile,
                                       region=self.region,
                                       numberOfProcessors=self.numberOfProcessors)


        try:
            num_reads_per_bin = np.concatenate(imap_res, axis=0)
            return num_reads_per_bin

        except ValueError:
            if self.bedFile:
                sys.exit('\nNo coverage values could be computed.\n\n'
                     'Please check that the chromosome names in the BED file are found on the bam files.\n\n'
                     'The valid chromosome names are:\n{}'.format(chrNames))
            else:
                sys.exit('\nNo coverage values could be computed.\n\nCheck that all bam files are valid and '
                     'contain mapped reads.')



    def count_reads_in_region(self, chrom, start, end, bed_regions_list=None):

        """Counts the reads in each bam file at each 'stepSize' position
        within the interval (start, end) for a window or bin of size binLength.

        The stepSize controls the distance between bins. For example,
        a step size of 20 and a bin size of 20 will create bins next to
        each other. If the step size is smaller than the bin size the
        bins will overlap.

        If a list of bedRegions is given, then the number of reads
        that overlaps with each region is counted.

        Parameters
        ----------
        chrom : str
            Chrom name
        start : int
            start coordinate
        end : int
            end coordinate
        bed_regions_list: list
            List of tuples of the form (chrom, start, end)
            corresponding to bed regions to be processed.
            If not bed file was passed to the object constructor
            then this list is empty.

        Returns
        -------
        numpy array
            The result is a numpy array that as rows each bin
            and as columns each bam file.


        Examples
        --------
        Initialize some useful values

        >>> test = Tester()
        >>> c = CountReadsPerBin([test.bamFile1, test.bamFile2], 25, 0, 0,
        ... stepSize=50)

        The transpose is used to get better looking numbers. The first line
        corresponds to the number of reads per bin in the first bamfile.

        >>> np.transpose(c.count_reads_in_region(test.chrom, 0, 200))
        array([[ 0.,  0.,  1.,  1.],
               [ 0.,  1.,  1.,  2.]])

        """

        if start > end:
            raise NameError("start %d bigger that end %d" % (start, end))

        if self.stepSize is None:
            raise ValueError("stepSize is not set!")
        # array to keep the read counts for the regions
        subNum_reads_per_bin = []

        rows = 0
        startTime = time.time()

        bamHandlers = [bamHandler.openBam(bam) for bam in self.bamFilesList]

        regionsToConsider = []
        if bed_regions_list is not None:
            for chrom, start, end in bed_regions_list:
                regionsToConsider.append((chrom, start, end, end - start))
        else:
            for i in xrange(start, end, self.stepSize):
                if i + self.binLength > end:
                    break
                regionsToConsider.append((chrom, i, i + self.binLength, self.binLength))

        for chrom, start, end, region_length in regionsToConsider:
            coverage_array = []
            for bam in bamHandlers:
                coverage_array.append(
                    self.get_coverage_of_region(bam, chrom, start, end, region_length)[0])

            subNum_reads_per_bin.extend(coverage_array)
            rows += 1

        if debug:
            endTime = time.time()
            print "%s countReadsInRegions_worker: processing %d " \
                  "(%.1f per sec) @ %s:%s-%s" % \
                  (multiprocessing.current_process().name,
                   rows, rows / (endTime - startTime), chrom, start, end )

        return np.array(subNum_reads_per_bin).reshape(rows, len(self.bamFilesList))


    def get_coverage_of_region(self, bamHandle, chrom, start, end, tileSize,
                            fragmentFromRead_func=None):

        """
        Returns a numpy array that corresponds to the number of reads
        that overlap with each tile.

        >>> test = Tester()
        >>> import pysam
        >>> c = CountReadsPerBin([], stepSize=1,
        ... defaultFragmentLength=300, extendPairedEnds=True)

        For this case the reads are length 36. The number of overlapping
        read fragments is 4 and 5 for the positions tested.

        >>> c.get_coverage_of_region(pysam.AlignmentFile(test.bamFile_PE), 'chr2',
        ... 5000833, 5000835, 1)
        array([ 4.,  5.])

        In the following example a paired read is extended to the fragment length which is 100
        The first mate starts at 5000000 and the second at 5000064. Each mate is
        extended to the fragment length *independently*
        At position 500090-500100 one fragment  of length 100 overlap, and after position 5000101
        there should be zero reads.

        >>> c.zerosToNans = True
        >>> c.get_coverage_of_region(pysam.AlignmentFile(test.bamFile_PE), 'chr2', 5000090, 5000110, 10)
        array([  1.,  nan])

        In the following  case the reads length is 50. Reads are not extended.

        >>> c.extendPairedEnds=False
        >>> c.defaultFragmentLength=1
        >>> c.get_coverage_of_region(pysam.AlignmentFile(test.bamFile2), '3R', 148, 154, 2)
        array([ 1.,  2.,  2.])


        """
        if not fragmentFromRead_func:
            fragmentFromRead_func = self.get_fragment_from_read
        length = end - start
        assert tileSize > 0, "bin length has to be an integer greater than zero. Current value {}".format(tileSize)
        if length % tileSize > 0:
            newLength = length - (length % tileSize)
            if debug:
                print  "length of region ({}) is not a multiple of " \
                       "tileSize {}\nThe region is being chopped to length " \
                       "{} bp".format(length, tileSize, newLength)

        vectorLength = length / tileSize
        coverage = np.zeros(vectorLength, dtype='float64')

        startTime = time.time()
        # caching seems faster. TODO: profile the function
        c = 0
        if chrom in bamHandle.references:
            # r.flag & 4 == 0 is to skip unmapped reads
            reads = [r for r in bamHandle.fetch(chrom, start, end)
                     if r.flag & 4 == 0]
        else:
            raise NameError("chromosome {} not found in bam file".format(chrom))

        prev_start_pos = None  # to store the start positions
        # of previous processed read pair
        for read in reads:
            if self.minMappingQuality and read.mapq < self.minMappingQuality:
                continue

            # filter reads based on SAM flag
            if self.samFlag_include and read.flag & self.samFlag_include == 0:
                continue
            if self.samFlag_exclude and read.flag & self.samFlag_exclude != 0:
                continue

            # get rid of duplicate reads that have same position on each of the
            # pairs
            if self.ignoreDuplicates and prev_start_pos \
                    and prev_start_pos == (read.reference_start, read.pnext, read.is_reverse):
                continue

            # since reads can be split (e.g. RNA-seq reads) each part of the
            # read that maps is called a position block.
            try:
                position_blocks = fragmentFromRead_func(read)
            except TypeError:
                # the get_fragment_from_read functions returns None in some cases.
                # Those cases are to be skipped, hence the continue line.
                continue

            for fragmentStart, fragmentEnd in position_blocks:
                fragmentLength = fragmentEnd - fragmentStart
                if fragmentLength == 0:
                    continue
                # skip reads that are not in the region being
                # evaluated.
                if fragmentEnd < start or fragmentStart >= end:
                    continue

                vectorStart = max((fragmentStart - start) / tileSize, 0)
                vectorEnd = min(np.ceil(float(fragmentEnd - start) / tileSize).astype('int'),
                                vectorLength)

                assert vectorEnd > vectorStart, "Error, vector end < than vector start"

                coverage[vectorStart:vectorEnd] += 1

            prev_start_pos = (read.reference_start, read.pnext, read.is_reverse)
            c += 1

        if debug:
            endTime = time.time()
            print "%s,  processing %s (%.1f per sec) reads @ %s:%s-%s" % (
                multiprocessing.current_process().name, c, c / (endTime - startTime), chrom, start, end)

        # change zeros to NAN
        if self.zerosToNans:
            coverage[coverage == 0] = np.nan

        return coverage

    def getReadLength(self, read):
        return len(read)


    def get_fragment_from_read(self, read):
        """Get read start and end position of a read.
        If given, the reads are extended as follows:
        If reads are paired end, each read mate is extended to match
        the fragment length, otherwise, a default fragment length
        is used. If reads are split (give by the CIGAR string) then
        the multiple positions of the read are returned.
        When reads are extended the cigar information is
        skipped.

        Parameters
        ----------
        read: pysam object.

        The following values are defined (for forward reads)::


                 |--          -- read.tlen --              --|
                 |-- read.alen --|
            -----|===============>------------<==============|----
                 |               |            |
            read.reference_start
                        read.reference_end  read.pnext

              and for reverse reads


                 |--             -- read.tlen --           --|
                                             |-- read.alen --|
            -----|===============>-----------<===============|----
                 |                           |               |
              read.pnext           read.reference_start  read.reference_end

        this is a sketch of a pair-end reads

        The function returns the fragment start and end, either
        using the paired end information (if available) or
        extending the read in the appropriate direction if this
        is single-end.

        Parameters
        ----------
        read : pysam read object


        Returns
        -------
        list of tuples
            [(fragment start, fragment end)]


        >>> test = Tester()
        >>> c = CountReadsPerBin([], 1, 1, 200)

        >>> c.extendPairedEnds = True
        >>> c.get_fragment_from_read(test.getRead("paired-forward"))
        [(5000000, 5000100)]
        >>> c.get_fragment_from_read(test.getRead("paired-reverse"))
        [(5000000, 5000100)]
        >>> c.get_fragment_from_read(test.getRead("single-forward"))
        [(5001491, 5001691)]
        >>> c.get_fragment_from_read(test.getRead("single-reverse"))
        [(5001536, 5001736)]
        >>> c.defaultFragmentLength = 0
        >>> c.get_fragment_from_read(test.getRead("single-forward"))
        [(5001491, 5001527)]
        >>> c.defaultFragmentLength = 0
        >>> c.extendPairedEnds = False
        >>> c.get_fragment_from_read(test.getRead("paired-forward"))
        [(5000000, 5000036)]

        Tests for read centering.

        >>> c.extendPairedEnds = True
        >>> c.center_read = True
        >>> c.defaultFragmentLength = 200
        >>> c.get_fragment_from_read(test.getRead("paired-forward"))
        [(5000032, 5000068)]
        >>> c.get_fragment_from_read(test.getRead("single-reverse"))
        [(5001618, 5001654)]
        """
        # if no extension is needed, use pysam get_blocks
        # to identify start and end reference positions.
        # get_blocks return a list of start and end positions
        # based on the CIGAR if skipped regions are found.
        # E.g for a cigar of 40M260N22M
        # get blocks return two elements for the first 40 matches
        # and the for the last 22 matches.
        if read.is_paired and self.extendPairedEnds == False:
            return read.get_blocks()
        if not read.is_paired and self.defaultFragmentLength <= read.alen:
            return read.get_blocks()

        if self.extendPairedEnds == True and read.is_paired \
                and abs(read.tlen) < self.maxPairedFragmentLength \
                and abs(read.tlen) > 0:
            if read.is_reverse:
                fragmentStart = read.pnext
                fragmentEnd = read.reference_end
            else:
                fragmentStart = read.reference_start
                # the end of the fragment is defined as
                # the start of the forward read plus the insert length
                fragmentEnd = read.reference_start + read.tlen

        else:
            if read.is_reverse:
                fragmentStart = read.reference_end - self.defaultFragmentLength
                fragmentEnd = read.reference_end
            else:
                fragmentStart = read.reference_start
                fragmentEnd = read.reference_start + self.defaultFragmentLength

        if self.center_read == True:
            fragmentCenter = fragmentEnd - (fragmentEnd - fragmentStart) / 2
            fragmentStart = fragmentCenter - read.alen / 2
            fragmentEnd = fragmentStart + read.alen

        return [(fragmentStart, fragmentEnd)]

    def getSmoothRange(self, tileIndex, tileSize, smoothRange, maxPosition):

        """
        Given a tile index position and a tile size (length), return the a new indices
        over a larger range, called the smoothRange.
        This region is centered in the tileIndex  an spans on both sizes
        to cover the smoothRange. The smoothRange is trimmed in case it is less
        than zero or greater than  maxPosition ::


             ---------------|==================|------------------
                        tileStart
                   |--------------------------------------|
                   |    <--      smoothRange     -->      |
                   |
             tileStart - (smoothRange-tileSize)/2

        Test for a smooth range that spans 3 tiles.

        Examples
        --------

        >>> c = CountReadsPerBin([], 1, 1, 1, 0)
        >>> c.getSmoothRange(5, 1, 3, 10)
        (4, 7)

        Test smooth range truncated on start.

        >>> c.getSmoothRange(0, 10, 30, 200)
        (0, 2)

        Test smooth range truncated on start.

        >>> c.getSmoothRange(1, 10, 30, 4)
        (0, 3)

        Test smooth range truncated on end.

        >>> c.getSmoothRange(5, 1, 3, 5)
        (4, 5)

        Test smooth range not multiple of tileSize.

        >>> c.getSmoothRange(5, 10, 24, 10)
        (4, 6)
        """
        smoothTiles = int(smoothRange / tileSize)
        if smoothTiles == 1:
            return (tileIndex, tileIndex + 1)

        smoothTilesSide = float(smoothTiles - 1) / 2
        smoothTilesLeft = int(np.ceil(smoothTilesSide))
        smoothTilesRight = int(np.floor(smoothTilesSide)) + 1

        indexStart = max(tileIndex - smoothTilesLeft, 0)
        indexEnd = min(maxPosition, tileIndex + smoothTilesRight)
        return (indexStart, indexEnd)


def remove_row_of_zeros(matrix):
    # remove rows containing all zeros or all nans
    _mat = np.nan_to_num(matrix)
    to_keep = _mat.sum(1) != 0
    return matrix[to_keep, :]

class Tester(object):
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
        self.bamFile1 = self.root + "testA.bam"
        self.bamFile2 = self.root + "testB.bam"
        self.bamFile_PE = self.root + "test_paired2.bam"
        self.chrom = '3R'
        global debug
        debug = 0

    def getRead(self, readType):
        """ prepare arguments for test
        """
        bam = bamHandler.openBam(self.bamFile_PE)
        if readType == 'paired-reverse':
            read = [x for x in bam.fetch('chr2', 5000081, 5000082)][0]
        elif readType == 'single-forward':
            read = [x for x in bam.fetch('chr2', 5001491, 5001492)][0]
        elif readType == 'single-reverse':
            read = [x for x in bam.fetch('chr2', 5001700, 5001701)][0]
        else:  # by default a forward paired read is returned
            read = [x for x in bam.fetch('chr2', 5000027, 5000028)][0]
        return read


