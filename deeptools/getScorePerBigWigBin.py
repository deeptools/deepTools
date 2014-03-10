import numpy as np
import time
import multiprocessing

# deepTools packages
import mapReduce
from bx.bbi.bigwig_file import BigWigFile

debug = 0


def countReadsInRegions_wrapper(args):
    return countReadsInRegions_worker(*args)


def countReadsInRegions_worker(chrom, start, end, bigwigFilesList,
                               bedRegions=None):
    """ returns the average score in each bigwig file at each
    'stepSize' position
    within the interval start, end for a 'binLength' window.
    Because the idea is to get counts for window positions at
    different positions for sampling the bins are equally spaced
    between each other and are  not one directly next *after* the other.

    If a list of bedRegions is given, then the number of reads
    that overlaps with each region is counted.

    """

    assert start < end, "start {} bigger that end {}".format(start, end)

    # array to keep the read counts for the regions
    subNum_reads_per_bin = []

    rows = 0

    bigWigHandlers = [BigWig(file=open(bw, 'r')) for bw in bigWigFilesList]
                      
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
        for bigWig in bigWigHandlers:
            score = bigWig.query(chrom, start, end, 1)[0]
            avgReadsArray.append(score['mean'])
        sub_score_per_bin.extend(avgReadsArray)
        rows += 1

    # the output is a matrix having as many rows as the variable 'row'
    # and as many columns as bigwig files. The rows correspond to
    # each of the regions processed by the worker.
    # np.array([[score1_1, score1_2],
    #           [score2_1, score2_2]]
    return np.array(sub_score_per_bin).reshape(rows, len(bigWigFilesList))


def getScorePerBin(bigwigFilesList, binLength, numberOfSamples,
                   numberOfProcessors=1,
                   verbose=False, region=None,
                   bedFile=None,
                   chrsToSkip=[]):

    r"""
    This function visits a number of sites and returs a matrix containing read
    counts. Each row to one sampled site and each column correspond to each of
    the bigwig files.

    """

    # Try to determine an optimal fraction of the genome (chunkSize)
    # that is sent to workers for analysis. If too short, too much time
    # is spend loading the files
    # if too long, some processors end up free.
    # the following values are empirical

    bigWigFilesHandlers = [BigWigFile(file=open(x, 'r'))
                           for x in bigWigFilesList]

    chromSizes = getCommonChrNames(bamFilesHandlers, verbose=verbose)

    # skip chromosome in the list. This is usually for the
    # X chromosome which may have either one copy  in a male sample
    # or a mixture of male/female and is unreliable.
    # Also the skip may contain heterochromatic regions and
    # mitochondrial DNA
    if len(chrsToSkip):
        chromSizes = [ x for x in chromSizes if x[0] not in chrsToSkip ]

    chrNames, chrLengths = zip(*chromSizes)

    genomeSize = sum(chrLengths)
    max_mapped = max([x.mapped for x in  bigwigFilesList])

    reads_per_bp = float(max_mapped) / genomeSize
#    chunkSize =  int(100 / ( reads_per_bp  * len(bamFilesList)) )

    stepSize = max(int( float(genomeSize) / numberOfSamples ), 1 )

    chunkSize =  int (stepSize * 1e3 / ( reads_per_bp  * len(bigwigFilesList)) )
    [bam_h.close() for bam_h in bigwigFilesList]

    if verbose:
        print "step size is {}".format(stepSize)

    if region:
        # in case a region is used, append the tilesize
        region += ":{}".format(binLength)

    imap_res = mapReduce.mapReduce((bigwigFilesList, stepSize, binLength),
                                    countReadsInRegions_wrapper,
                                    chromSizes,
                                    genomeChunkLength=chunkSize,
                                    bedFile=bedFile,
                                    region=region,
                                    numberOfProcessors=numberOfProcessors)

    score_per_bin = np.concatenate(imap_res, axis=0)
    return score_per_bin


def getReadLength(read):
    return len(read)


def getCoverageOfRegion(bamHandle, chrom, start, end, tileSize,
                        defaultFragmentLength, extendPairedEnds=True, 
                        zerosToNans=True, maxPairedFragmentLength=None,
                        minMappingQuality=None, ignoreDuplicates=False,
                        fragmentFromRead_func = getFragmentFromRead):
    """
    Returns a numpy array that corresponds to the number of reads 
    that overlap with each tile.

    >>> test = Tester()
    >>> import pysam

    For this case the reads are length 36. For the positions given
    the number of overlapping read fragments is 4 and 5
    >>> getCoverageOfRegion(pysam.Samfile(test.bamFile_PE), 'chr2', 5000833, 5000835, 1, 0, False)
    array([ 4.,  5.])

    In the following example a paired read is extended to the fragment length wich is 100
    The first mate starts at 5000000 and the second at 5000064. Each mate is
    extended to the fragment length *independently*
    At position 500090-500100 one fragment  of length 100 overlap, and after position 5000101  
    there should be zero reads
    >>> getCoverageOfRegion(pysam.Samfile(test.bamFile_PE), 'chr2', 5000090, 5000110, 10, 0, True)
    array([  1.,  nan])

    In the following  case the reads length is 50
    >>> getCoverageOfRegion(pysam.Samfile(test.bamFile2), '3R', 148, 154, 2, 0, False)
    array([ 1.,  2.,  2.])

    Test ignore duplicates
    >>> getCoverageOfRegion(pysam.Samfile(test.bamFile2), '3R', 0, 200, 50, 0, False, ignoreDuplicates=True)
    array([ nan,   1.,   1.,   1.])

    Test long regions
    >>> getCoverageOfRegion(pysam.Samfile(test.bamFile2), '3R', 0, 200, 200, 0, False)
    array([ 4.])
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
                 if r.flag & 4 == 0]
    else:
        raise NameError("chromosome {} not found in bam file".format(chrom))

    prev_start_pos = None # to store the start positions
                          # of previous processed read pair

    for read in reads:
        if minMappingQuality and read.mapq < minMappingQuality:
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
                                      maxPairedFragmentLength=maxPairedFragmentLength)

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
    than cero or greater than  maxPosition


     ---------------|==================|------------------
                tileStart
           |--------------------------------------|
           |    <--      smoothRange     -->      |
           |        
     tileStart - (smoothRange-tileSize)/2

    Test for a smooth range that spans 3 tiles
    >>> getSmoothRange(5, 1, 3, 10)
    (4, 7)
    
    Test smooth range truncated on start
    >>> getSmoothRange(0, 10, 30, 200)
    (0, 2)

    Test smooth range truncated on start
    >>> getSmoothRange(1, 10, 30, 4)
    (0, 3)

    Test smooth range truncated on end
    >>> getSmoothRange(5, 1, 3, 5)
    (4, 5)

    Test smooth range not multiple of tileSize
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


