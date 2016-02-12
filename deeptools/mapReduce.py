import multiprocessing

debug = 0


def mapReduce(staticArgs, func, chromSize,
              genomeChunkLength=None,
              region=None,
              bedFile=None,
              blackListFileName=None,
              numberOfProcessors=4,
              verbose=False,
              self_=None):
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

    :param chromSize: A list of duples containing the chromosome
                      name and its length
    :param region: The format is chr:start:end:tileSize (see function
                   getUserRegion)
    :param staticArgs: tuple of arguments that are sent to the given 'func'

    :param func: function to call. The function is called using the
                 following parameters (chrom, start, end, staticArgs)
    :param bedFile: Is a bed file is given, the args to the func to be
                    called are extended to include a list of bed
                    defined regions.
    :param blackListFileName: A list of regions to exclude from all computations.
                              Note that this has genomeChunkLength resolution...
    :param self_: In case mapreduce should make a call to an object
                the self variable has to be passed.
    """

    if not genomeChunkLength:
        genomeChunkLength = 1e5

    if verbose:
        print "genome partition size for multiprocessing: {}".format(
            genomeChunkLength)

    region_start = 0
    region_end = None

    # if a region is set, that means that the task should be only cover
    # the given genomic position

    if region:
        chromSize, region_start, region_end, genomeChunkLength = getUserRegion(chromSize, region)
        if verbose:
            print "chrom size: {}, region start: {}, region end: {}, " \
                  "genome chunk length sent to each procesor: {}".format(chromSize, region_start, region_end, genomeChunkLength)

    if bedFile:
        bed_interval_tree = BED_to_interval_tree(bedFile)
        # modify chromSize such that it only contains
        # chromosomes that are in the bed file
        if not region:
            chromSize = [x for x in chromSize if x[0] in bed_interval_tree.keys()]
            if not len(chromSize):
                exit("*ERROR*\nChromosome names in bed file do not match the chromosome names files to process")

        else:
            bed_in_region = bed_interval_tree[chromSize[0][0]].find(region_start, region_end)
            if chromSize[0][0] not in bed_interval_tree.keys() or len(bed_in_region) == 0:
                exit("*ERROR*\nThe specified region {} does not contain any of the "
                     "intervals in the bed file".format(region))
            # the user region has to be extended, otherwise if a bed interval ends after the region is not counted
            # since the user region has been already added to the chromSize, the chrom size is extended instead.
            if bed_in_region[-1].end > chromSize[0][1]:
                chromSize[0] = (chromSize[0][0], bed_in_region[-1].end)

    if blackListFileName:
        blackList = BED_to_interval_tree(open(blackListFileName, "r"))

    TASKS = []
    # iterate over all chromosomes
    for chrom, size in chromSize:
        # the start is zero unless a specific region is defined
        start = 0 if region_start == 0 else region_start
        for startPos in xrange(start, size, genomeChunkLength):
            endPos = min(size, startPos + genomeChunkLength)

            # Reject a chunk if it overlaps
            if blackListFileName:
                regions = blSubtract(blackList, chrom, [startPos, endPos])
            else:
                regions = [[startPos, endPos]]

            for reg in regions:
                if self_ is not None:
                    argsList = [self_]
                else:
                    argsList = []

                argsList.extend([chrom, reg[0], reg[1]])
                # add to argument list the static list received the the function
                argsList.extend(staticArgs)

                # if a bed file is given, append to the TASK list,
                # a list of bed regions that overlap with the
                # current genomeChunk.
                if bedFile:
                    # this method to get the bedFile regions may seem
                    # cumbersome but I (fidel) think is better to
                    # balance the load between multiple processors.
                    # This method first partitions the genome into smaller
                    # chunks and then, for each chunk, the list of
                    # regions overlapping the chunk interval is added.
                    # This is preferable to sending each worker a
                    # single region because of the overhead of initiating
                    # the data.
                    bed_regions_list = []
                    for bed_region in bed_interval_tree[chrom].find(reg[0], reg[1]):
                        # start + 1 is used to avoid regions that may overlap
                        # with two genomeChunks to be counted twice. Such region
                        # is only added for the genomeChunk that contains the start
                        # of the bed region.
                        if bed_region.start < endPos < bed_region.end:
                            continue
                        bed_regions_list.append([chrom, bed_region.start, bed_region.end])
                    if len(bed_regions_list) == 0:
                        continue
                    # add to argument list, the position of the bed regions to use
                    argsList.append(bed_regions_list)

                TASKS.append(tuple(argsList))

    if len(TASKS) > 1 and numberOfProcessors > 1:
        if verbose:
            print ("using {} processors for {} "
                   "number of tasks".format(numberOfProcessors,
                                            len(TASKS)))

        pool = multiprocessing.Pool(numberOfProcessors)
        res = pool.map_async(func, TASKS).get(9999999)
    else:
        res = map(func, TASKS)

    return res


def getUserRegion(chrom_sizes, region_string, max_chunk_size=1e6):
    r"""
    Verifies if a given region argument, given by the user
    is valid. The format of the region_string is chrom:start:end:tileSize
    where start, end and tileSize are optional.

    :param chrom_sizes: dictionary of chromosome/scaffold size. Key=chromosome name
    :param region_string: a string of the form chr:start:end
    :param max_chunk_size: upper limit for the chunk size
    :return: tuple chrom_size for the region start, region end, chunk size

    #>>> data = getUserRegion({'chr2': 1000}, "chr1:10:10")
    #Traceback (most recent call last):
    #    ...
    #NameError: Unknown chromosome: chr1
    #Known chromosomes are: ['chr2']

    If the region end is biger than the chromosome size, this
    value is used instead
    >>> getUserRegion({'chr2': 1000}, "chr2:10:1001")
    ([('chr2', 1000)], 10, 1000, 990)

    Test chunk and regions size reduction to match tile size
    >>> getUserRegion({'chr2': 200000}, "chr2:10:123344:3")
    ([('chr2', 123344)], 9, 123345, 123336)
    """
    region = region_string.split(":")
    chrom = region[0]
    chrom_sizes = dict(chrom_sizes)

    try:
        chrom_sizes[chrom]
    except KeyError:
        raise NameError("Unknown chromosome: %s\nKnown "
                        "chromosomes are: %s " % (chrom, chrom_sizes.keys()))
    try:
        region_start = int(region[1])
    except IndexError:
        region_start = 0
    try:
        region_end = int(region[2]) if int(region[2]) <= chrom_sizes[chrom] \
            else chrom_sizes[chrom]
    except IndexError:
        region_end = chrom_sizes[chrom]
    if region_start > region_end or region_start < 0:
        raise NameError("{} not valid. The format is chrom:start:end. "
                        "Without comas, dashes or dots. ".format(region_string))
    try:
        tilesize = int(region[3])
    except IndexError:
        tilesize = None

    chrom_sizes = [(chrom, region_end)]

    # if tilesize is given, make region_start and region_end
    # multiple of tileSize
    if tilesize:
        region_start -= region_start % tilesize
        region_end += tilesize - (region_end % tilesize)

    chunk_size = int(region_end - region_start)
    if chunk_size > max_chunk_size:
        chunk_size = max_chunk_size
        if tilesize and tilesize < chunk_size:
            chunk_size -= chunk_size % tilesize

    return chrom_sizes, region_start, region_end, int(chunk_size)


def BED_to_interval_tree(BED_file):
    """
    Creates an index of intervals, using an interval tree, for each BED entry

    :param BED_file: file handler of a BED file

    :return interval tree
    """
    from bx.intervals.intersection import IntervalTree, Interval

    bed_interval_tree = {}
    for line in BED_file:
        if line[0] == "#":
            continue
        fields = line.strip().split()
        chrom, start_bed, end_bed, = fields[0], int(fields[1]), int(fields[2])

        if chrom not in bed_interval_tree.keys():
            bed_interval_tree[chrom] = IntervalTree()

        # skip if a region overlaps with a region already seen
        """
        if len(bed_interval_tree[chrom].find(start_bed, start_bed + 1)) > 0:
            continue
        """
        bed_interval_tree[chrom].add_interval(Interval(start_bed, end_bed))

    return bed_interval_tree


def blOverlap(t, chrom, chunk):
    """
    Test for an overlap between an IntervalTree and a given genomic chunk.

    This attempts to account for differences in chromosome naming.

    Returns the overlaps
    """

    if t is None:
        return []

    if chrom not in t.keys():
        if chrom.startswith("chr"):
            chrom = chrom[3:]
        elif chrom == "MT":
            chrom = "chrM"
        else:
            chrom = "chr" + chrom

        if chrom not in t.keys():
            return []

    return t[chrom].find(chunk[0], chunk[1])


def blSubtract(t, chrom, chunk):
    """
    If a genomic region overlaps with a blacklisted region, then subtract that region out

    returns a list of lists
    """

    if t is None:
        return [chunk]

    overlaps = blOverlap(t, chrom, chunk)
    if len(overlaps) > 0:
        output = []
        for o in overlaps:
            if chunk[1] <= chunk[0]:
                break
            if chunk[0] < o.start:
                output.append([chunk[0], o.start])
            chunk[0] = o.end
        if chunk[0] < chunk[1]:
            output.append([chunk[0], chunk[1]])
    else:
        output = [chunk]

    return output
