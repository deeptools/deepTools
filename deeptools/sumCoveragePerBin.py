import numpy as np
import multiprocessing
import time

from deeptools import countReadsPerBin
from deeptools.utilities import getTLen
from deeptoolsintervals import GTF


class SumCoveragePerBin(countReadsPerBin.CountReadsPerBin):
    r"""This is an extension of CountReadsPerBin for use with plotFingerprint.
    There, we need to sum the per-base coverage.
    """
    def get_coverage_of_region(self, bamHandle, chrom, regions,
                               fragmentFromRead_func=None):
        """
        Returns a numpy array that corresponds to the number of reads
        that overlap with each tile.

        >>> test = Tester()
        >>> import pysam
        >>> c = SumCoveragePerBin([], stepSize=1, extendReads=300)

        For this case the reads are length 36. The number of overlapping
        read fragments is 4 and 5 for the positions tested. Note that reads are
        NOT extended, due to there being a 0 length input list of BAM files!

        >>> c.get_coverage_of_region(pysam.AlignmentFile(test.bamFile_PE), 'chr2',
        ... [(5000833, 5000834), (5000834, 5000835)])
        array([4., 5.])

        In the following  case the reads length is 50. Reads are not extended.

        >>> c.extendReads=False
        >>> c.get_coverage_of_region(pysam.AlignmentFile(test.bamFile2), '3R', [(148, 150), (150, 152), (152, 154)])
        array([2., 4., 4.])


        """
        if not fragmentFromRead_func:
            fragmentFromRead_func = self.get_fragment_from_read
        nbins = len(regions)
        if len(regions[0]) == 3:
            nbins = 0
            for reg in regions:
                nbins += (reg[1] - reg[0]) // reg[2]
        coverages = np.zeros(nbins, dtype='float64')

        if self.defaultFragmentLength == 'read length':
            extension = 0
        else:
            extension = self.maxPairedFragmentLength

        blackList = None
        if self.blackListFileName is not None:
            blackList = GTF(self.blackListFileName)

        vector_start = 0
        for idx, reg in enumerate(regions):
            if len(reg) == 3:
                tileSize = int(reg[2])
                nRegBins = (reg[1] - reg[0]) // tileSize
            else:
                nRegBins = 1
                tileSize = int(reg[1] - reg[0])

            # Blacklisted regions have a coverage of 0
            if blackList and blackList.findOverlaps(chrom, reg[0], reg[1]):
                continue
            regStart = int(max(0, reg[0] - extension))
            regEnd = reg[1] + int(extension)

            # If alignments are extended and there's a blacklist, ensure that no
            # reads originating in a blacklist are fetched
            if blackList and reg[0] > 0 and extension > 0:
                o = blackList.findOverlaps(chrom, regStart, reg[0])
                if o is not None and len(o) > 0:
                    regStart = o[-1][1]
                o = blackList.findOverlaps(chrom, reg[1], regEnd)
                if o is not None and len(o) > 0:
                    regEnd = o[0][0]

            start_time = time.time()
            # caching seems faster. TODO: profile the function
            c = 0
            try:
                # BAM input
                if chrom not in bamHandle.references:
                    raise NameError("chromosome {} not found in bam file".format(chrom))
            except:
                # bigWig input, as used by plotFingerprint
                if bamHandle.chroms(chrom):
                    _ = np.array(bamHandle.stats(chrom, regStart, regEnd, type="mean", nBins=nRegBins), dtype=np.float)
                    _[np.isnan(_)] = 0.0
                    _ = _ * tileSize
                    coverages += _
                    continue
                else:
                    raise NameError("chromosome {} not found in bigWig file with chroms {}".format(chrom, bamHandle.chroms()))

            prev_pos = set()
            lpos = None
            # of previous processed read pair
            for read in bamHandle.fetch(chrom, regStart, regEnd):
                if read.is_unmapped:
                    continue
                if self.minMappingQuality and read.mapq < self.minMappingQuality:
                    continue

                # filter reads based on SAM flag
                if self.samFlag_include and read.flag & self.samFlag_include != self.samFlag_include:
                    continue
                if self.samFlag_exclude and read.flag & self.samFlag_exclude != 0:
                    continue

                # Fragment lengths
                tLen = getTLen(read)
                if self.minFragmentLength > 0 and tLen < self.minFragmentLength:
                    continue
                if self.maxFragmentLength > 0 and tLen > self.maxFragmentLength:
                    continue

                # get rid of duplicate reads that have same position on each of the
                # pairs
                if self.ignoreDuplicates:
                    # Assuming more or less concordant reads, use the fragment bounds, otherwise the start positions
                    if tLen >= 0:
                        s = read.pos
                        e = s + tLen
                    else:
                        s = read.pnext
                        e = s - tLen
                    if read.reference_id != read.next_reference_id:
                        e = read.pnext
                    if lpos is not None and lpos == read.reference_start \
                            and (s, e, read.next_reference_id, read.is_reverse) in prev_pos:
                        continue
                    if lpos != read.reference_start:
                        prev_pos.clear()
                    lpos = read.reference_start
                    prev_pos.add((s, e, read.next_reference_id, read.is_reverse))

                # since reads can be split (e.g. RNA-seq reads) each part of the
                # read that maps is called a position block.
                try:
                    position_blocks = fragmentFromRead_func(read)
                except TypeError:
                    # the get_fragment_from_read functions returns None in some cases.
                    # Those cases are to be skipped, hence the continue line.
                    continue

                last_eIdx = None
                for fragmentStart, fragmentEnd in position_blocks:
                    if fragmentEnd is None or fragmentStart is None:
                        continue
                    fragmentLength = fragmentEnd - fragmentStart
                    if fragmentLength == 0:
                        continue
                    # skip reads that are not in the region being
                    # evaluated.
                    if fragmentEnd <= reg[0] or fragmentStart >= reg[1]:
                        continue

                    if fragmentStart < reg[0]:
                        fragmentStart = reg[0]
                    if fragmentEnd > reg[0] + len(coverages) * tileSize:
                        fragmentEnd = reg[0] + len(coverages) * tileSize

                    sIdx = vector_start + max((fragmentStart - reg[0]) // tileSize, 0)
                    eIdx = vector_start + min(np.ceil(float(fragmentEnd - reg[0]) / tileSize).astype('int'), nRegBins)
                    if eIdx >= len(coverages):
                        eIdx = len(coverages) - 1
                    if last_eIdx is not None:
                        sIdx = max(last_eIdx, sIdx)
                        if sIdx >= eIdx:
                            continue

                    # First bin
                    if fragmentEnd < reg[0] + (sIdx + 1) * tileSize:
                        _ = fragmentEnd - fragmentStart
                    else:
                        _ = reg[0] + (sIdx + 1) * tileSize - fragmentStart
                    if _ > tileSize:
                        _ = tileSize
                    coverages[sIdx] += _
                    _ = sIdx + 1
                    while _ < eIdx:
                        coverages[_] += tileSize
                        _ += 1
                    while eIdx - sIdx >= nRegBins:
                        eIdx -= 1
                    if eIdx > sIdx:
                        _ = fragmentEnd - (reg[0] + eIdx * tileSize)
                        if _ > tileSize:
                            _ = tileSize
                        elif _ < 0:
                            _ = 0
                        coverages[eIdx] += _
                    last_eIdx = eIdx

                c += 1

            if self.verbose:
                endTime = time.time()
                print("%s,  processing %s (%.1f per sec) reads @ %s:%s-%s" % (
                    multiprocessing.current_process().name, c, c / (endTime - start_time), chrom, reg[0], reg[1]))

            vector_start += nRegBins

        # change zeros to NAN
        if self.zerosToNans:
            coverages[coverages == 0] = np.nan

        return coverages


class Tester(object):

    def __init__(self):
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
        import os
        self.root = os.path.dirname(os.path.abspath(__file__)) + "/test/test_data/"
        self.bamFile1 = self.root + "testA.bam"
        self.bamFile2 = self.root + "testB.bam"
        self.bamFile_PE = self.root + "test_paired2.bam"
        self.chrom = '3R'
