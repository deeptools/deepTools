import os

import pyBigWig

debug = 0


def bedGraphToBigWig(chromSizes, bedGraphFiles, bigWigPath):
    """
    Takes a sorted list of bedgraph files and write them to a single bigWig file using pyBigWig.
    The order of bedGraphFiles must match that of chromSizes!
    """
    bw = pyBigWig.open(bigWigPath, "w")
    assert bw is not None
    bw.addHeader(chromSizes, maxZooms=10)
    lastChrom = None
    starts = []
    ends = []
    vals = []
    for bg in bedGraphFiles:
        if bg is not None:
            f = open(bg)
            for line in f:
                interval = line.split()
                # Buffer up to a million entries
                if interval[0] != lastChrom or len(starts) == 1000000:
                    if lastChrom is not None:
                        bw.addEntries([lastChrom] * len(starts), starts, ends=ends, values=vals)
                    lastChrom = interval[0]
                    starts = [int(interval[1])]
                    ends = [int(interval[2])]
                    vals = [float(interval[3])]
                else:
                    starts.append(int(interval[1]))
                    ends.append(int(interval[2]))
                    vals.append(float(interval[3]))
            f.close()
            os.remove(bg)
    if len(starts) > 0:
        bw.addEntries([lastChrom] * len(starts), starts, ends=ends, values=vals)
    bw.close()


def getGenomeChunkLength(bamHandles, tile_size, mappedList):
    """
    Tries to estimate the length of the genome sent to the workers
    based on the density of reads per bam file and the number
    of bam files.

    The chunk length should be a multiple of the tileSize

    """

    genomeLength = sum(bamHandles[0].lengths)

    max_reads_per_bp = max([float(x) / genomeLength for x in mappedList])

    # 2e6 is an empirical estimate
    genomeChunkLength = int(min(5e6, int(2e6 / (max_reads_per_bp * len(bamHandles)))))

    genomeChunkLength -= genomeChunkLength % tile_size
    return genomeChunkLength
