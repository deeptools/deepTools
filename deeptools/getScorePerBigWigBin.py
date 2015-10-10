import pyBigWig
import numpy as np
import os
import sys
import warnings

# deepTools packages
import deeptools.mapReduce as mapReduce
import deeptools.config as cfg
#debug = 0


def countReadsInRegions_wrapper(args):
    # Using arguments unpacking!
    return countFragmentsInRegions_worker(*args)


def countFragmentsInRegions_worker(chrom, start, end,
                               bigWigFiles,
                               stepSize, binLength,     #staticArgs
                               bedRegions=None):
    """ returns the average score in each bigwig file at each 'stepSize'
    position within the interval start, end for a 'binLength' window.
    Because the idea is to get counts for window positions at
    different positions for sampling the bins are equally spaced
    and *not adjacent*.

    If a list of bedRegions is given, then the number of reads
    that overlaps with each region is counted.

    Test dataset with two samples covering 200 bp.
    >>> test = Tester()

    Fragment coverage.
    >>> np.transpose(countFragmentsInRegions_worker(test.chrom, 0, 200,
    ... [test.bwFile1, test.bwFile2], 50, 25))
    array([[ 1.,  1.,  2.,  2.],
           [ 1.,  1.,  1.,  3.]])

    >>> np.transpose(countFragmentsInRegions_worker(test.chrom, 0, 200,
    ... [test.bwFile1, test.bwFile2], 200, 200))
    array([[ 1.5],
           [ 1.5]])

    BED regions:
    >>> bedRegions = [(test.chrom, 45, 55), (test.chrom, 95, 105),
    ... (test.chrom, 145, 155)]
    >>> np.transpose(countFragmentsInRegions_worker(test.chrom, 0, 200,
    ... [test.bwFile1, test.bwFile2], 200, 200, bedRegions=bedRegions))
    array([[ 1. ,  1.5,  2. ],
           [ 1. ,  1. ,  2. ]])
    """

    assert start < end, "start {} bigger that end {}".format(start, end)

    # array to keep the scores for the regions
    sub_score_per_bin = []

    rows = 0

    regionsToConsider = []
    if bedRegions:
        for chrom, start, end in bedRegions:
            regionsToConsider.append((chrom, start, end, end - start))
    else:
        for i in xrange(start, end, stepSize):
            if (i + binLength) > end:
                regionsToConsider.append((chrom, i, end, end-i)) #last bin (may be smaller)
            else:
                regionsToConsider.append((chrom, i, i + binLength, binLength))

    #print "\nNumber of regions: {}".format(len(regionsToConsider))

    warnings.simplefilter("default")
    i = 0
    num_warnings = 0
    for chrom, start, end, binLength in regionsToConsider:
        avgReadsArray = []
        i += 1
        for idx, bwh in enumerate(bigWigFiles):
            score = bwh.stats(chrom, start, end) #The default is "mean" and 1 bin

            if score is None or np.isnan(score):
                if num_warnings < 10:
                    sys.stderr.write("NaN found in {} at {}:{:,}-{:,}\n".format(bigwigFilesList[idx], chrom, start, end))
                    num_warnings += 1
                score = np.nan
            avgReadsArray.append(score)     #mean of fragment coverage for region
        #print "{} Region: {}:{:,}-{:,} {}  {} {}".format(i, chrom, start, end, binLength, avgReadsArray[0], avgReadsArray[1])

        sub_score_per_bin.extend(avgReadsArray)
        rows += 1
    warnings.resetwarnings()

    # the output is a matrix having as many rows as the variable 'row'
    # and as many columns as bigwig files. The rows correspond to
    # each of the regions processed by the worker.
    # np.array([[score1_1, score1_2],
    #           [score2_1, score2_2]]

    return np.array(sub_score_per_bin).reshape(rows, len(bigwigFilesList))


def getChromSizes(bigWigFiles):
    """
    Get chromosome sizes from bigWig files.

    Test dataset with two samples covering 200 bp.
    >>> test = Tester()

    Chromosome name(s) and size(s).
    >>> getChromSizes([test.bwFile1, test.bwFile2])
    [('3R', 200)]
    """

    chromNamesAndSize = {}
    for bw in bigWigFiles :
        if(bw is None or !bw) :
            return None

        for k,v in bw.chroms() :
            if(k not in chromNamesAndSize) :
                chromNamesAndSize[k] = v
            else if(chromNamesAndSize[k] != v) :
                print "\nWARNING\n" \
                "Chromosome {} length reported in the bigwig files differ.\n\n" \
                "The smaller of the two will be used.".format(k, chomNamesAndSize[k], v)
                if(chromNamesAndSize[k] >= v) :
                    chromNamesAndSize[k] = v
    # get the list of common chromosome names and sizes
    chromSizes = sorted([(k, v) for k, v in chromNamesAndSize.iteritems() ])
    return chromSizes


#This is a terribly inefficient way to get this.
def getNumberOfFragmentsPerRegionFromBigWig(bw, chromSizes):
    """
    Get the number of all mapped fragments per region in all chromosomes
    from a bigWig. Utilizing bx-python.

    Test dataset with two samples covering 200 bp.
    >>> test = Tester()

    Get number of fragments in sample.
    >>> getNumberOfFragmentsPerRegionFromBigWig(test.bwFile1, [('3R', 200)])
    3.0
    >>> getNumberOfFragmentsPerRegionFromBigWig(test.bwFile2, [('3R', 200)])
    4.0
    """
    mapped = 0
    for cname, csize in chromSizes:
        regions = bw.intervals(cname, 0, csize) # region = bwh.get(chrom_name, start, end)
        for region in regions:
            mapped += region[2]
    return mapped


def getScorePerBin(bigWigFiles, binLength,
                   numberOfProcessors=1,
                   verbose=False, region=None,
                   bedFile=None,
                   stepSize=None,
                   chrsToSkip=[]):
    """
    This function returns a matrix containing scores (median) for the coverage
    of fragments within a region. Each row corresponds to a sampled region.
    Likewise, each column corresponds to a bigwig file.

    Test dataset with two samples covering 200 bp.
    >>> test = Tester()

    >>> np.transpose(getScorePerBin([test.bwFile1, test.bwFile2], 50, 5,))
    array([[ 1.,  1.,  2.,  2.],
           [ 1.,  1.,  1.,  3.]])

    """

    # Try to determine an optimal fraction of the genome (chunkSize)
    # that is sent to workers for analysis. If too short, too much time
    # is spent loading the files
    # if too long, some processors end up free.
    # the following is a heuristic

    # get list of common chromosome names and sizes
    chromSizes = getChromSizes(bigWigFiles)

    # skip chromosome in the list. This is usually for the
    # X chromosome which may have either one copy  in a male sample
    # or a mixture of male/female and is unreliable.
    # Also the skip may contain heterochromatic regions and
    # mitochondrial DNA
    if len(chrsToSkip): chromSizes = [ x for x in chromSizes if x[0] not in chrsToSkip ]

    chrNames, chrLengths = zip(*chromSizes)
    genomeSize = sum(chrLengths)
    if stepSize is None:
        stepSize = binLength    #for adjacent bins
    chunkSize = int(stepSize * 500 / len(bigWigFiles))
    if verbose:
        print "step size is {}".format(stepSize)

    if region:
        # in case a region is used, append the tilesize
        region += ":{}".format(binLength)
    # mapReduce( (staticArgs), func, chromSize, etc. )
    imap_res = mapReduce.mapReduce((bigWigFiles, stepSize, binLength),
                                    countReadsInRegions_wrapper,
                                    chromSizes,
                                    genomeChunkLength=chunkSize,
                                    bedFile=bedFile,
                                    region=region,
                                    numberOfProcessors=numberOfProcessors)

    score_per_bin = np.concatenate(imap_res, axis=0)
    return score_per_bin


class Tester():
    def __init__( self ):
        """
        The distribution of reads (and fragments) in the two bigWig
        files is as follows.

        They cover 200 bp::

              0                              100                           200
              |------------------------------------------------------------|
            A                                ==============>- - - - - - - -
              - - - - - - - - - - - - - - - - - - - - - - - <==============


            B - - - - - - - - <==============               ==============>
                                             ==============>- - - - - - - -
                                                            ==============>
        """
        self.root = "./test/test_data/"
        self.bwFile1  = self.root + "testA.bw"
        self.bwFile2  = self.root + "testB.bw"
        self.bwFile_PE  = self.root + "test_paired2.bw"
        self.chrom = '3R'
        #global debug
        #debug = 0
