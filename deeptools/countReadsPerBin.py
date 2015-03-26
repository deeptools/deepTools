import numpy as np
import time
import multiprocessing

# deepTools packages
from utilities import getCommonChrNames
import bamHandler
import mapReduce

debug = 0


def countReadsInRegions_wrapper(args):
    """
    Passes the arguments to countReadsInRegions_worker.
    This is a step to make the code more readable, given
    the constrains from the multiprocessing module.
    """

    return countReadsInRegions_worker(*args)


def countReadsInRegions_worker(chrom, start, end, bamFilesList,
                               stepSize, binLength, defaultFragmentLength,
                               skipZeros=False,
                               extendPairedEnds=True,
                               minMappingQuality=None,
                               ignoreDuplicates=False,
                               samFlag=None,
                               bedRegions=None
                               ):
    """Counts the reads in each bam file at each 'stepSize' position
    within the interval (start, end) for a window or bin of size binLength.
    Because the idea is to get counts for window/bin positions at

    The stepSize controls the distance between bins. For example,
    a step size of 20 and a bin size of 20 will create bins next to
    each other. if the step size is smaller than the bin size the
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
    bamFileList : list
        List of name of indexed bam files.
    stepSize : int
        the positions for which the coverage is computed are defined as follows:
        ``range(start, end, stepSize)``. Thus, a stepSize of 1, will compute
        the coverage at each base pair. If the stepSize is equal to the
        binLength then the coverage is computed for consecutive bins. If seepSize is
        smaller than the binLength, then teh bins will overlap.
    binLength : int
        length of the window/bin
    defaultFragmentLength : in
        see :meth:`deepTools.countReadsPerBin.getFragmentFromRead` method

    Returns
    -------
    numpy array
        The result is a numpy array that as rows each bin
        and as columns each bam file.


    Examples
    --------
    Initialize some useful values

    >>> test = Tester()

    The transpose is used to get better looking numbers. The first line
    corresponds to the number of reads per bin in the first bamfile.

    >>> np.transpose(countReadsInRegions_worker(test.chrom, 0, 200,
    ... [test.bamFile1, test.bamFile2], 50, 25, 0))
    array([[ 0.,  0.,  1.,  1.],
           [ 0.,  1.,  1.,  2.]])

    When skipZeros is set to true, those cases in which *all* of the
    bamfiles have zero counts for a certain bin are ignored.

    >>> np.transpose(countReadsInRegions_worker(test.chrom, 0, 200,
    ... [test.bamFile1, test.bamFile2], 50, 25, 0, skipZeros=True))
    array([[ 0.,  1.,  1.],
           [ 1.,  1.,  2.]])


    """

    if start > end:
        raise NameError("start %d bigger that end %d" % (start, end))

    # array to keep the read counts for the regions
    subNum_reads_per_bin = []

    rows = 0
    startTime = time.time()
    extendPairedEnds = True
    zerosToNans = False

    bamHandlers = [bamHandler.openBam(bam) for bam in bamFilesList]

    regionsToConsider = []
    if bedRegions:
        for chrom, start, end in bedRegions:
            regionsToConsider.append((chrom, start, end, end - start))
    else:
        for i in xrange(start, end, stepSize):
            if i + binLength > end:
                break
            regionsToConsider.append((chrom, i, i + binLength, binLength))

    for chrom, start, end, binLength in regionsToConsider:
        avgReadsArray = []
        for bam in bamHandlers:
            avgReadsArray.append(
                getCoverageOfRegion(bam,
                                    chrom, start, end,
                                    binLength,
                                    defaultFragmentLength,
                                    extendPairedEnds,
                                    zerosToNans,
                                    minMappingQuality=minMappingQuality,
                                    ignoreDuplicates=ignoreDuplicates,
                                    samFlag=samFlag
                                    )[0])
        # skip if any of the bam files returns a NaN
        if np.isnan(sum(avgReadsArray)):
            continue

        if skipZeros and sum(avgReadsArray) == 0:
            continue
        subNum_reads_per_bin.extend(avgReadsArray)
        rows += 1

    if debug:
        endTime = time.time()
        print "%s countReadsInRegions_worker: processing %d " \
            "(%.1f per sec) @ %s:%s-%s"  % \
            (multiprocessing.current_process().name,
             rows, rows / (endTime - startTime), chrom, start, end )

    return np.array(subNum_reads_per_bin).reshape(rows, len(bamFilesList))


def getNumReadsPerBin(bamFilesList, binLength, numberOfSamples,
                      defaultFragmentLength, numberOfProcessors=1,
                      skipZeros=True, verbose=False, region=None,
                      bedFile=None, extendPairedEnds=True,
                      minMappingQuality=None,
                      ignoreDuplicates=False,
                      chrsToSkip=[],
                      stepSize=None,
                      samFlag=None):

    r"""
    This function collects read counts (coverage) from several bam files and returns
    an numpy array with the results. This function does not explicitly do the
    coverage computation, instead divides the work into smaller chunks that are
    sent to individual processors.

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

    skipZeros : bool
        Default is True. This option decides if regions having zero coverage in all bam files
        should be skipped or kept.

    verbose : bool
        Output messages. Default: False

    region : str
        Region to limit the computation in the form chrom:start:end.

    bedFile : str
        Name of a bed file containing the regions for wich to compute the coverage. This option
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

    samFlag : int
        If given, only reads having such flag are considered. For example, to get only
        reads that are the first mates a samFlag of 64 could be used. Similarly, the
        samFlag can be used to select only reads mapping on the forward (or reverse) strand
        or to get only properly paired reads.

    Returns
    -------
    numpy array

        Each row correspond to each bin/bed region and each column correspond to each of
        the bamFiles. If ``skipZeros`` is used, then the result may have less rows
        than expected


    Examples
    --------

    The test data contains reads for 200 bp.

    >>> test = Tester()

    The transpose function is used to get a nicer looking output.
    The first line corresponds to the number of reads per bin in bam file 1

    >>> np.transpose(getNumReadsPerBin([test.bamFile1, test.bamFile2],
    ... 50, 4, 0, skipZeros=True))
    array([[ 0.,  1.,  1.],
           [ 1.,  1.,  2.]])
    """

    # Try to determine an optimal fraction of the genome (chunkSize) that is sent to 
    # workers for analysis. If too short, too much time is spend loading the files
    # if too long, some processors end up free.
    # the following values are empirical

    bamFilesHandlers = [ bamHandler.openBam(x) for x in bamFilesList ]
    chromSizes = getCommonChrNames(bamFilesHandlers, verbose=verbose)

    # skip chromosome in the list. This is usually for the
    # X chromosome which may have either one copy  in a male sample
    # or a mixture of male/female and is unreliable.
    # Also the skip may contain heterochromatic regions and
    # mitochondrial DNA
    if len(chrsToSkip): chromSizes = [ x for x in chromSizes if x[0] not in chrsToSkip ]

    chrNames, chrLengths = zip(*chromSizes)

    genomeSize = sum(chrLengths)
    max_mapped = max( [ x.mapped for x in  bamFilesHandlers ] )

    reads_per_bp = float(max_mapped) / genomeSize
#    chunkSize =  int(100 / ( reads_per_bp  * len(bamFilesList)) )

    if stepSize is None:
        stepSize = max(int( float(genomeSize) / numberOfSamples ), 1 )

    chunkSize =  int (stepSize * 1e3 / ( reads_per_bp  * len(bamFilesHandlers)) )
    [ bam_h.close() for bam_h in bamFilesHandlers]

    if verbose:
        print "step size is {}".format(stepSize)

    if region:
        # in case a region is used, append the tilesize
        region += ":{}".format(binLength)

    imap_res = mapReduce.mapReduce( (bamFilesList, stepSize, binLength,
                                     defaultFragmentLength, skipZeros,
                                     extendPairedEnds, minMappingQuality,
                                     ignoreDuplicates, samFlag),
                                    countReadsInRegions_wrapper,
                                    chromSizes,
                                    genomeChunkLength=chunkSize,
                                    bedFile=bedFile,
                                    region=region,
                                    numberOfProcessors=numberOfProcessors)

    try:
        num_reads_per_bin = np.concatenate(imap_res, axis=0)
    except ValueError:
        if bedFile:
            exit('\nNo coverage values could be computed.\n\n'
                 'Please check that the chromosome names in the BED file are found on the bam files.\n\n'
                 'The valid chromosome names are:\n{}'.format(chrNames))
        else:
            exit('\nNo coverage values could be computed.\n\nCheck that all bam files are valid and '
                 'contain mapped reads.')

    return num_reads_per_bin


def getReadLength(read):
    return len(read)


def getFragmentFromRead(read, defaultFragmentLength, extendPairedEnds=True, 
                        maxPairedFragmentLength=None, centerRead=False):
    """
    The read has to be pysam object.

    The following values are defined (for forward reads)::


             |--          -- read.tlen --              --|
             |-- read.alen --|
        -----|===============>------------<==============|----
             |               |            |
          read.pos      read.aend      read.pnext


          and for reverse reads


             |--             -- read.tlen --           --|
                                         |-- read.alen --|
        -----|===============>-----------<===============|----
             |                           |               |
          read.pnext                   read.pos      read.aend

    this is a sketch of a pair-end reads

    The function returns the fragment start and end, either
    using the paired end information (if available) or
    extending the read in the appropriate direction if this
    is single-end.

    >>> test = Tester()
    >>> 
    >>> getFragmentFromRead(test.getRead("paired-forward"), 200)
    (5000000, 5000100)
    >>> getFragmentFromRead(test.getRead("paired-reverse"), 200)
    (5000000, 5000100L)
    >>> getFragmentFromRead(test.getRead("single-forward"), 200)
    (5001491, 5001691)
    >>> getFragmentFromRead(test.getRead("single-reverse"), 200)
    (5001536L, 5001736L)
    >>> getFragmentFromRead(test.getRead("single-forward"), 20)
    (5001491, 5001527L)
    >>> getFragmentFromRead(test.getRead("paired-forward"), 0, False)
    (5000000, 5000036L)

    Tests for read centering.

    >>> getFragmentFromRead(test.getRead("paired-forward"), 200,
    ... True, centerRead=True)
    (5000032, 5000068)
    >>> getFragmentFromRead(test.getRead("single-reverse"), 200,
    ... centerRead=True)
    (5001618L, 5001654L)
    """
    # convert reads to fragments

    fragmentStart = fragmentEnd = None
    # this option indicates that the paired ends correspond
    # to the fragment ends
    # condition read.tlen < maxPairedFragmentLength is added to avoid read pairs
    # that spand thousands of base pairs
    if not maxPairedFragmentLength:
        maxPairedFragmentLength = 2*defaultFragmentLength if defaultFragmentLength > 0 else 1000;

    if extendPairedEnds == True and read.is_paired  \
            and abs(read.tlen) < maxPairedFragmentLength \
            and abs(read.tlen) > 0:
        if read.is_reverse:
            fragmentStart = read.pnext
            fragmentEnd   = read.aend
        else:
            fragmentStart = read.pos
            # the end of the fragment is defined as 
            # the start of the forward read plus the insert length
            fragmentEnd   = read.pos + read.tlen
    else:
        if defaultFragmentLength <= read.aend - read.pos:
            fragmentStart = read.pos
            fragmentEnd   = read.aend
        else:
            if read.is_reverse:
                fragmentStart = read.aend - defaultFragmentLength
                fragmentEnd   = read.aend
            else:
                fragmentStart = read.pos
                fragmentEnd   = read.pos + defaultFragmentLength

    if centerRead == True:
        fragmentCenter = fragmentEnd - (fragmentEnd - fragmentStart)/2
        fragmentStart = fragmentCenter - read.alen/2
        fragmentEnd   = fragmentStart + read.alen

    return (fragmentStart, fragmentEnd)


def getCoverageOfRegion(bamHandle, chrom, start, end, tileSize,
                        defaultFragmentLength, extendPairedEnds=True, 
                        zerosToNans=True, maxPairedFragmentLength=None,
                        minMappingQuality=None, ignoreDuplicates=False,
                        fragmentFromRead_func=getFragmentFromRead,
                        centerRead=False, samFlag=None):
    """
    Returns a numpy array that corresponds to the number of reads 
    that overlap with each tile.

    >>> test = Tester()
    >>> import pysam

    For this case the reads are length 36. The number of overlapping
    read fragments is 4 and 5 for the positions tested.

    >>> getCoverageOfRegion(pysam.Samfile(test.bamFile_PE), 'chr2',
    ... 5000833, 5000835, 1, 0, False)
    array([ 4.,  5.])

    In the following example a paired read is extended to the fragment length which is 100
    The first mate starts at 5000000 and the second at 5000064. Each mate is
    extended to the fragment length *independently*
    At position 500090-500100 one fragment  of length 100 overlap, and after position 5000101  
    there should be zero reads.

    >>> getCoverageOfRegion(pysam.Samfile(test.bamFile_PE), 'chr2', 5000090, 5000110, 10, 0, True)
    array([  1.,  nan])

    In the following  case the reads length is 50.

    >>> getCoverageOfRegion(pysam.Samfile(test.bamFile2), '3R', 148, 154, 2, 0, False)
    array([ 1.,  2.,  2.])

    Test ignore duplicates.

    >>> getCoverageOfRegion(pysam.Samfile(test.bamFile2), '3R', 0, 200, 50, 0,
    ... False, ignoreDuplicates=True)
    array([ nan,   1.,   1.,   1.])

    Test long regions.

    >>> getCoverageOfRegion(pysam.Samfile(test.bamFile2), '3R', 0, 200, 200, 0, False)
    array([ 4.])

    Test sam flag with value = 64 which means only first mate.

    >>> getCoverageOfRegion(pysam.Samfile(test.bamFile_PE), 'chr2', 5000833,
    ... 5000835, 1, 0, False, samFlag=64)
    array([ nan,   1.])

    """
    if not fragmentFromRead_func:
        fragmentFromRead_func = getFragmentFromRead
    length = end - start
    if length % tileSize > 0:
        newLength = length - (length % tileSize)
        if debug:
            print  "length of region ({}) is not a multiple of "\
                "tileSize {}\nThe region is being chopped to length "\
                "{} bp".format(length, tileSize, newLength)

    vectorLength = length / tileSize
    coverage = np.zeros(vectorLength, dtype='float64')

    startTime = time.time()
    # caching seems faster. TODO: profile the function
    c = 0
    if chrom in bamHandle.references:
        # r.flag & 4 == 0 is to skip unmapped reads
        reads = [r for r in bamHandle.fetch(chrom, start, end)
                 if r.flag & 4 == 0 ]
    else:
        raise NameError("chromosome {} not found in bam file".format(chrom))

    prev_start_pos = None # to store the start positions
                          # of previous processed read pair

    for read in reads:
        if minMappingQuality and read.mapq < minMappingQuality:
            continue

        # filter reads based on SAM flag
        if samFlag and read.flag & samFlag == 0:
            continue

        # get rid of duplicate reads that have same position on each of the
        # pairs
        if ignoreDuplicates and prev_start_pos  \
                and prev_start_pos == (read.pos, read.pnext, read.is_reverse):
            continue

        try:
            fragmentStart, fragmentEnd = \
                fragmentFromRead_func(read, defaultFragmentLength,
                                      extendPairedEnds,
                                      maxPairedFragmentLength=maxPairedFragmentLength,
                                      centerRead=centerRead)

            fragmentLength = fragmentEnd - fragmentStart
            if fragmentLength == 0:
                fragmentLength = defaultFragmentLength

            vectorStart = max( (fragmentStart - start)/tileSize, 0)
            vectorEnd   = min( np.ceil(float(fragmentEnd   - start)/tileSize).astype('int'), 
                               vectorLength)

            coverage[vectorStart:vectorEnd] +=  1
        except TypeError:
            # the getFragmentFromRead functions returns None in some cases.
            # Those cases are to be skiped, hence the continue line.
            continue

        prev_start_pos = (read.pos, read.pnext, read.is_reverse)
        c += 1

    if debug:
        endTime = time.time()
        print "%s,  processing %s (%.1f per sec) reads @ %s:%s-%s" % (multiprocessing.current_process().name, c, c / (endTime - startTime) ,chrom, start, end)

    # change zeros to NAN
    if zerosToNans:
        coverage[coverage == 0] = np.nan

    return coverage 

def getSmoothRange(tileIndex, tileSize, smoothRange, maxPosition):
    
    """
    Given a tile index position and a tile size (length), return the a new indices
    over a larger range, called the smoothRange.
    This region is centered in the tileIndex  an spans on both sizes
    to cover the smoothRange. The smoothRange is trimed in case it is less
    than zero or greater than  maxPosition ::


         ---------------|==================|------------------
                    tileStart
               |--------------------------------------|
               |    <--      smoothRange     -->      |
               |
         tileStart - (smoothRange-tileSize)/2

    Test for a smooth range that spans 3 tiles.

    >>> getSmoothRange(5, 1, 3, 10)
    (4, 7)
    
    Test smooth range truncated on start.

    >>> getSmoothRange(0, 10, 30, 200)
    (0, 2)

    Test smooth range truncated on start.

    >>> getSmoothRange(1, 10, 30, 4)
    (0, 3)

    Test smooth range truncated on end.

    >>> getSmoothRange(5, 1, 3, 5)
    (4, 5)

    Test smooth range not multiple of tileSize.

    >>> getSmoothRange(5, 10, 24, 10)
    (4, 6)
    """
    smoothTiles = int(smoothRange/tileSize)
    if smoothTiles == 1:
        return (tileIndex, tileIndex + 1)

    smoothTilesSide = float(smoothTiles - 1) / 2
    smoothTilesLeft = int(np.ceil(smoothTilesSide))
    smoothTilesRight = int(np.floor(smoothTilesSide)) + 1

    indexStart = max( tileIndex - smoothTilesLeft, 0 )
    indexEnd   = min( maxPosition, tileIndex + smoothTilesRight )
    return (indexStart, indexEnd)

class Tester():
    def __init__( self ):
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

    def getRead( self, readType ):
        """ prepare arguments for test
        """
        bam = bamHandler.openBam(self.bamFile_PE)
        if readType == 'paired-reverse':
            read = [x for x in bam.fetch('chr2', 5000081,5000082 )][0]
        elif readType == 'single-forward':
            read = [x for x in bam.fetch('chr2', 5001491, 5001492 )][0]
        elif readType == 'single-reverse':
            read = [x for x in bam.fetch('chr2', 5001700, 5001701 )][0]
        else: # by default a forward paired read is returned
            read = [x for x in bam.fetch('chr2', 5000027,5000028 )][0]
        return read


