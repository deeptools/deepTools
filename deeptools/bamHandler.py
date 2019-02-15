import sys
import pysam
from deeptools.mapReduce import mapReduce


def countReadsInInterval(args):
    chrom, start, end, fname, toEOF = args

    bam = openBam(fname)
    mapped = 0
    unmapped = 0
    for b in bam.fetch(chrom, start, end):
        if chrom == "*":
            unmapped += 1
            continue
        if b.pos < start:
            continue
        if not b.is_unmapped:
            mapped += 1
        else:
            unmapped += 1
    return mapped, unmapped, chrom


def getMappingStats(bam, nThreads):
    """
    This is used for CRAM files, since idxstats() and .mapped/.unmapped are meaningless

    This requires pysam > 0.13.0
    """
    header = [(x, y) for x, y in zip(bam.references, bam.lengths)]
    res = mapReduce([bam.filename, False], countReadsInInterval, header, numberOfProcessors=nThreads)

    mapped = sum([x[0] for x in res])
    unmapped = sum([x[1] for x in res])
    stats = {x[0]: [0, 0] for x in header}
    for r in res:
        stats[r[2]][0] += r[0]
        stats[r[2]][1] += r[1]

    # We need to count the number of unmapped reads as well
    unmapped += bam.count("*")

    return mapped, unmapped, stats


def openBam(bamFile, returnStats=False, nThreads=1, minimalDecoding=True):
    """
    A wrapper for opening BAM/CRAM files.

    bamFile: str
        A BAM/CRAM file name

    returnStats: bool
        Return a tuple of (file_handle, nMappedReads, nUnmappedReads, statsDict).
        These additional values are needed by some downstream functions, since one
        can't use file_handle.mapped on CRAM files (or idxstats())

    nThreads: int
        If returnStats is True, number of threads to use for computing statistics

    minimalDecoding: Bool
        For CRAM files, don't decode the read name, sequence, qual, or auxiliary tag fields (these aren't used by most functions).

    Returns either the file handle or a tuple as described in returnStats
    """
    format_options = ["required_fields=0x1FF"]
    if sys.version_info.major >= 3:
        format_options = [b"required_fields=0x1FF"]
    if not minimalDecoding:
        format_options = None
    try:
        bam = pysam.Samfile(bamFile, 'rb', format_options=format_options)
    except IOError:
        sys.exit("The file '{}' does not exist".format(bamFile))
    except:
        sys.exit("The file '{}' does not have BAM or CRAM format ".format(bamFile))

    try:
        assert(bam.check_index() is not False)
    except:
        sys.exit("'{}' does not appear to have an index. You MUST index the file first!".format(bamFile))

    if bam.is_cram and returnStats:
        mapped, unmapped, stats = getMappingStats(bam, nThreads)
    elif bam.is_bam:
        mapped = bam.mapped
        unmapped = bam.unmapped

        # Make the dictionary to hold the stats
        if returnStats:
            stats = {chrom.contig: [chrom.mapped, chrom.unmapped] for chrom in bam.get_index_statistics()}

    if bam.is_bam or (bam.is_cram and returnStats):
        if mapped == 0:
            sys.stderr.write("WARNING! '{}' does not have any mapped reads. Please "
                             "check that the file is properly indexed and "
                             "that it contains mapped reads.\n".format(bamFile))

    if returnStats:
        return bam, mapped, unmapped, stats
    else:
        return bam
