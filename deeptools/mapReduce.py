import multiprocessing

debug = 0
def mapReduce(staticArgs, func, chromSize,
              genomeChunkLength=None, 
              region=None,
              numberOfProcessors=4,
              verbose=False):

    """
    Split the genome into parts that are sent to workers using a defined
    number of procesors. Results are collected and returned.

    For each genomic region the given 'func' is called using
    the following parameters:

     chrom, start, end, staticArgs

    The *arg* are static, *pickable* variables that need to be sent
    to workers.

    The genome chunk length corresponds to a fraction of the genome, in bp,
    that is send to each of the workers for processing.

    Depending on the type of process a larger or shorter regions may be 
    preferred

    :param chromSize: A list of duples containing the chromome name and its length
    :param region: The format is chr:start:end
    """


    if not genomeChunkLength:
        genomeChunkLength = 1e5

    if verbose:
        print "genome partition size for multiprocessing: %s" % (genomeChunkLength) 

    regionStart = 0

    # if a region is set, that means that the task should be only cover
    # the given genomic possition

    if region:
        chromSize, regionStart, regionEnd, genomeChunkLength = getUserRegion(chromSize, region)
        if verbose:
            print (chromSize, regionStart, regionEnd, genomeChunkLength)

    TASKS = []
    # iterate over all chromosomes
    for chrom, size in chromSize:
        # the start is zero unless a specific region is defined
        start = 0 if regionStart == 0 else regionStart
        for startPos in xrange(start, size, genomeChunkLength):
            endPos = min(size, startPos + genomeChunkLength)
            argsList = [chrom, startPos, endPos]
            # add to argument list the static list received the the function
            argsList.extend(staticArgs)
            TASKS.append( tuple(argsList) )

    if len(TASKS) > 1 and numberOfProcessors>1:
        if verbose:
            print ("using {} processors for {} "
                   "number of tasks".format(numberOfProcessors,
                                            len(TASKS)))

        pool = multiprocessing.Pool(numberOfProcessors)
        res = pool.map_async( func, TASKS ).get(9999999)
    else:
        res = map(func, TASKS)

    return res


def getUserRegion( chromSizes, regionString,
                   max_chunk_size=1e6):
    """
    Verifies if a given region argument, given by the user
    is valid. The format of the regionString is chrom:start:end:tileSize
    where start, end and tileSize are optional.
    >>> data = getUserRegion({'chr2': 1000}, "chr1:10:10")
    Traceback (most recent call last):
    NameError: Unkown chromosome: chr1
    Known chromosomes are: ['chr2'] 

    If the region end is biger than the chromosome size, this value is used instead
    >>> getUserRegion({'chr2': 1000}, "chr2:10:1001")
    ([('chr2', 1000)], 10, 1000, 990)

    Test chunk and regions size reduction to match tile size
    >>> getUserRegion({'chr2': 200000}, "chr2:10:123344:3")
    ([('chr2', 123344)], 9, 123345, 123336)
    """
    region = regionString.split(":")
    chrom = region[0]
    chromSizes = dict(chromSizes)

    try:
        chromSizes[chrom]
    except KeyError:
        raise NameError("Unkown chromosome: %s\nKnown chromosomes are: %s " % (chrom, chromSizes.keys()))
    try: 
        regionStart = int(region[1])         
    except IndexError:
        regionStart = 0
    try:
        regionEnd = int(region[2]) if int(region[2]) <= chromSizes[chrom] else chromSizes[chrom]
    except IndexError:
        regionEnd = chromSizes[chrom]
    if regionStart > regionEnd or regionStart < 0:
        raise NameError( "%s not valid. The format is chrom:start:end. Without comas, dashes or dots. " % (regionString) )
    try:
        tilesize = int(region[3])
    except IndexError:
        tilesize = None

    chromSizes = [ (chrom, regionEnd) ]

    # if tilesize is given, make regionStart and regionEnd
    # multiple of tileSize
    if tilesize:
        regionStart -= regionStart  % tilesize
        regionEnd += tilesize - (regionEnd  % tilesize)
    
#    chunkSize = int((regionEnd - regionStart) / ( 2 * multiprocessing.cpu_count() ))
    chunkSize = int(regionEnd - regionStart)
    if chunkSize > max_chunk_size: 
        chunkSize=max_chunk_size
        if tilesize and tilesize < chunkSize:
            chunkSize -= chunkSize % tilesize

    return (chromSizes, regionStart, regionEnd, int(chunkSize))
