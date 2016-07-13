import multiprocessing
from deeptoolsintervals import GTF
import random

debug = 0


def mapReduce(staticArgs, func, chromSize,
              genomeChunkLength=None,
              region=None,
              bedFile=None,
              blackListFileName=None,
              numberOfProcessors=4,
              verbose=False,
              includeLabels=False,
              keepExons=False,
              transcriptID="transcriptID",
              exonID="exonID",
              transcript_id_designator="transcript_id",
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
    :param includeLabels: Pass group and transcript labels into the calling
                          function. These are added to the static args
                          (groupLabel and transcriptName).

    If "includeLabels" is true, a tuple of (results, labels) is returned
    """

    if not genomeChunkLength:
        genomeChunkLength = 1e5
    genomeChunkLength = int(genomeChunkLength)

    if verbose:
        print("genome partition size for multiprocessing: {0}".format(
            genomeChunkLength))

    region_start = 0
    region_end = None

    # if a region is set, that means that the task should be only cover
    # the given genomic position

    if region:
        chromSize, region_start, region_end, genomeChunkLength = getUserRegion(chromSize, region)
        if verbose:
            print("chrom size: {0}, region start: {1}, region end: {2}, "
                  "genome chunk length sent to each procesor: {3}".format(chromSize, region_start, region_end, genomeChunkLength))

    if bedFile:
        defaultGroup = None
        if len(bedFile) == 1:
            defaultGroup = "genes"
        bed_interval_tree = GTF(bedFile, defaultGroup=defaultGroup, transcriptID=transcriptID, exonID=exonID, transcript_id_designator=transcript_id_designator, keepExons=keepExons)

    if blackListFileName:
        blackList = GTF(blackListFileName)

    TASKS = []
    # iterate over all chromosomes
    for chrom, size in chromSize:
        # the start is zero unless a specific region is defined
        start = 0 if region_start == 0 else region_start
        for startPos in range(start, size, genomeChunkLength):
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
                    # This effectively creates batches of intervals, which is
                    # generally more performant due to the added overhead of
                    # initializing additional workers.

                    # TODO, there's no point in including the chromosome
                    if includeLabels:
                        bed_regions_list = [[chrom, x[4], x[2], x[3], x[5], x[6]] for x in bed_interval_tree.findOverlaps(chrom, reg[0], reg[1], trimOverlap=True, numericGroups=True, includeStrand=True)]
                    else:
                        bed_regions_list = [[chrom, x[4], x[5], x[6]] for x in bed_interval_tree.findOverlaps(chrom, reg[0], reg[1], trimOverlap=True, includeStrand=True)]

                    if len(bed_regions_list) == 0:
                        continue
                    # add to argument list, the position of the bed regions to use
                    argsList.append(bed_regions_list)

                TASKS.append(tuple(argsList))

    if len(TASKS) > 1 and numberOfProcessors > 1:
        if verbose:
            print(("using {} processors for {} "
                   "number of tasks".format(numberOfProcessors,
                                            len(TASKS))))
        random.shuffle(TASKS)
        pool = multiprocessing.Pool(numberOfProcessors)
        res = pool.map_async(func, TASKS).get(9999999)
    else:
        res = list(map(func, TASKS))

    if includeLabels:
        if bedFile:
            return res, bed_interval_tree.labels
        else:
            return res, None
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

    Test chromosome name mismatch
    >>> getUserRegion({'2': 200000}, "chr2:10:123344:3")
    ([('2', 123344)], 9, 123345, 123336)
    >>> getUserRegion({'chrM': 200000}, "MT:10:123344:3")
    ([('chrM', 123344)], 9, 123345, 123336)
    """
    region = region_string.split(":")
    chrom = region[0]
    chrom_sizes = dict(chrom_sizes)

    if chrom not in list(chrom_sizes.keys()):
        if chrom == "MT":
            chromUse = "chrM"
        elif chrom == "chrM":
            chromUse = "MT"
        elif chrom[0:3] == "chr":
            chromUse = chrom[3:]
        else:
            chromUse = "chr" + chrom
        if chromUse not in list(chrom_sizes.keys()):
            raise NameError("Unknown chromosome: %s\nKnown "
                            "chromosomes are: %s " % (chrom, list(chrom_sizes.keys())))
        chrom = chromUse
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


def blSubtract(t, chrom, chunk):
    """
    If a genomic region overlaps with a blacklisted region, then subtract that region out

    returns a list of lists
    """

    if t is None:
        return [chunk]

    overlaps = t.findOverlaps(chrom, chunk[0], chunk[1])
    if overlaps is not None and len(overlaps) > 0:
        output = []
        for o in overlaps:
            if chunk[1] <= chunk[0]:
                break
            if chunk[0] < o[0]:
                output.append([chunk[0], o[0]])
            chunk[0] = o[1]
        if chunk[0] < chunk[1]:
            output.append([chunk[0], chunk[1]])
    else:
        output = [chunk]

    return output
