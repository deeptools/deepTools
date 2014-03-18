from bx.bbi.bigwig_file import BigWigFile
import numpy as np
import os
import sys
import warnings

# deepTools packages
import mapReduce

#debug = 0


def countReadsInRegions_wrapper(args):
    # Using arguments unpacking!
    return countFragmentsInRegions_worker(*args)


def countFragmentsInRegions_worker(chrom, start, end,
                               bigwigFilesList,
                               stepSize, binLength,     #staticArgs
                               skipZeros=False,
                               bedRegions=None):
    """ returns the average score in each bigwig file at each 'stepSize'
    position within the interval start, end for a 'binLength' window.
    Because the idea is to get counts for window positions at
    different positions for sampling the bins are equally spaced
    between each other and are not one directly next *after* the other.

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

    bigWigHandlers = [BigWigFile(open(bw, "rb")) for bw in bigwigFilesList]

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

    for chrom, start, end, binLength in regionsToConsider:
        avgReadsArray = []
        i += 1
        for bwh in bigWigHandlers:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                score = bwh.query(chrom, start, end, 1)[0]
            if np.isnan(score['mean']) or score is None:
                sys.stderr.write("{}  {} found at {}:{:,}-{:,}\n".format(i, score['mean'], chrom, start, end))
                score['mean'] = 0.0
            avgReadsArray.append(score['mean'])     #mean of fragment coverage for region
        #print "{} Region: {}:{:,}-{:,} {}  {} {}".format(i, chrom, start, end, binLength, avgReadsArray[0], avgReadsArray[1])

        if skipZeros and sum(avgReadsArray) == 0:
            continue
        sub_score_per_bin.extend(avgReadsArray)
        rows += 1
    warnings.resetwarnings()

    # the output is a matrix having as many rows as the variable 'row'
    # and as many columns as bigwig files. The rows correspond to
    # each of the regions processed by the worker.
    # np.array([[score1_1, score1_2],
    #           [score2_1, score2_2]]

    return np.array(sub_score_per_bin).reshape(rows, len(bigwigFilesList))


def getChromSizes(bigwigFilesList):
    """
    Get chromosome sizes from bigWig file by shell calling bigWigInfo
    (UCSC tools).

    Test dataset with two samples covering 200 bp.
    >>> test = Tester()

    Chromosome name(s) and size(s).
    >>> getChromSizes([test.bwFile1, test.bwFile2])
    [('3R', 200)]
    """
    #The following lines are - with one exception ("bigWigInfo") -
    #identical with the bw-reading part of deeptools/countReadsPerBin.py (FK)
    cCommon = []
    chromNamesAndSize = {}
    for bw in bigwigFilesList:
        inBlock = False
        for line in os.popen("{} -chroms {}".format("bigWigInfo", bw)).readlines():
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
                            bw[0], size, bw[1])
                        chromNamesAndSize[chromName] = min(
                            chromNamesAndSize[chromName], size)
                else:
                    chromNamesAndSize[chromName] = size
    # get the list of common chromosome names and sizes
    chromSizes = sorted([(k, v) for k, v in chromNamesAndSize.iteritems() \
                         if k in cCommon])
    return chromSizes


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
    bwh = BigWigFile(open(bw, "rb"))
    mapped = 0
    for cname, csize in chromSizes:
        regions = bwh.get(cname, 0, csize) # region = bwh.get(chrom_name, start, end)
        for region in regions:
            mapped += region[2]
    return mapped


def getScorePerBin(bigwigFilesList, binLength,
                   numberOfProcessors=1, skipZeros=True,
                   verbose=False, region=None,
                   bedFile=None,
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
    # is spend loading the files
    # if too long, some processors end up free.
    # the following values are empirical

    # get list of common chromosome names and sizes
    chromSizes = getChromSizes(bigwigFilesList)

    # skip chromosome in the list. This is usually for the
    # X chromosome which may have either one copy  in a male sample
    # or a mixture of male/female and is unreliable.
    # Also the skip may contain heterochromatic regions and
    # mitochondrial DNA
    if len(chrsToSkip): chromSizes = [ x for x in chromSizes if x[0] not in chrsToSkip ]

    chrNames, chrLengths = zip(*chromSizes)
    genomeSize = sum(chrLengths)
    max_mapped = max( map(lambda x: getNumberOfFragmentsPerRegionFromBigWig(x, chromSizes), bigwigFilesList) )
    reads_per_bp = float(max_mapped) / genomeSize
    stepSize = binLength    #for consecutive bins
    chunkSize = int(stepSize * 1e3 / ( reads_per_bp  * len(bigwigFilesList)) )

    if verbose:
        print "step size is {}".format(stepSize)

    if region:
        # in case a region is used, append the tilesize
        region += ":{}".format(binLength)

    # mapReduce( (staticArgs), func, chromSize, etc. )
    imap_res = mapReduce.mapReduce((bigwigFilesList, stepSize, binLength, skipZeros),
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

        They cover 200 bp.

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
