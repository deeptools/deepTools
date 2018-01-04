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
            print([b.query_name, b.pos])
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


def openBam(bamFile, returnStats=False, nThreads=1):
    try:
        bam = pysam.Samfile(bamFile, 'rb')
    except IOError:
        sys.exit("The file '{}' does not exist".format(bamFile))
    except:
        sys.exit("The file '{}' does not have BAM format ".format(bamFile))

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
            sys.exit("'{}' does not have any mapped reads. Please "
                     "check that the file is properly indexed and "
                     "that it contains mapped reads.".format(bamFile))

    if returnStats:
        return bam, mapped, unmapped, stats
    else:
        return bam
