import shutil
import os
import time
import sys
import multiprocessing
import numpy as np

# deepTools packages
import deeptools.utilities
from deeptools import bamHandler
from deeptools import mapReduce
from deeptoolsintervals import GTF

debug = 0
old_settings = np.seterr(all='ignore')


def countReadsInRegions_wrapper(args):
    """
    Passes the arguments to countReadsInRegions_worker.
    This is a step required given
    the constrains from the multiprocessing module.
    The args var, contains as first element the 'self' value
    from the countReadsPerBin object

    """
    return CountReadsPerBin.count_reads_in_region(*args)


class CountReadsPerBin(object):

    r"""Collects coverage over multiple bam files using multiprocessing

    This function collects read counts (coverage) from several bam files and returns
    an numpy array with the results. This class uses multiprocessing to compute the coverage.

    Parameters
    ----------
    bamFilesList : list
        List containing the names of indexed bam files. E.g. ['file1.bam', 'file2.bam']

    binLength : int
        Length of the window/bin. This value is overruled by ``bedFile`` if present.

    numberOfSamples : int
        Total number of samples. The genome is divided into ``numberOfSamples``, each
        with a window/bin length equal to ``binLength``. This value is overruled
        by ``stepSize`` in case such value is present and by ``bedFile`` in which
        case the number of samples and bins are defined in the bed file

    numberOfProcessors : int
        Number of processors to use. Default is 4

    verbose : bool
        Output messages. Default: False

    region : str
        Region to limit the computation in the form chrom:start:end.

    bedFile : file_handle
        File handle of a bed file containing the regions for which to compute the coverage. This option
        overrules ``binLength``, ``numberOfSamples`` and ``stepSize``.

    blackListFileName : str
        A string containing a BED file with blacklist regions.

    extendReads : bool, int

        Whether coverage should be computed for the extended read length (i.e. the region covered
        by the two mates or the regions expected to be covered by single-reads).
        If the value is 'int', then then this is interpreted as the fragment length to extend reads
        that are not paired. For Illumina reads, usual values are around 300.
        This value can be determined using the peak caller MACS2 or can be
        approximated by the fragment lengths computed when preparing the library for sequencing. If the value
        is of the variable is true and not value is given, the fragment size is sampled from the library but
        only if the library is paired-end. Default: False


    minMappingQuality : int
        Reads of a mapping quality less than the give value are not considered. Default: None

    ignoreDuplicates : bool
        Whether read duplicates (same start, end position. If paired-end, same start-end for mates) are
        to be excluded. Default: false

    chrToSkip: list
        List with names of chromosomes that do not want to be included in the coverage computation.
        This is useful to remove unwanted chromosomes (e.g. 'random' or 'Het').

    stepSize : int
        the positions for which the coverage is computed are defined as follows:
        ``range(start, end, stepSize)``. Thus, a stepSize of 1, will compute
        the coverage at each base pair. If the stepSize is equal to the
        binLength then the coverage is computed for consecutive bins. If seepSize is
        smaller than the binLength, then teh bins will overlap.

    center_read : bool
        Determines if reads should be centered with respect to the fragment length.

    samFlag_include : int
        Extracts only those reads having the SAM flag. For example, to get only
        reads that are the first mates a samFlag of 64 could be used. Similarly, the
        samFlag_include can be used to select only reads mapping on the reverse strand
        or to get only properly paired reads.

    samFlag_exclude : int
        Removes reads that match the SAM flag. For example to get all reads
        that map to the forward strand a samFlag_exlude 16 should be used. Which
        translates into exclude all reads that map to the reverse strand.

    zerosToNans : bool
        If true, zero values encountered are transformed to Nans. Default false.

    minFragmentLength : int
        If greater than 0, fragments below this size are excluded.

    maxFragmentLength : int
        If greater than 0, fragments above this size are excluded.

    out_file_for_raw_data : str
        File name to save the raw counts computed

    Returns
    -------
    numpy array

        Each row correspond to each bin/bed region and each column correspond to each of
        the bamFiles.


    Examples
    --------

    The test data contains reads for 200 bp.

    >>> test = Tester()

    The transpose function is used to get a nicer looking output.
    The first line corresponds to the number of reads per bin in bam file 1

    >>> c = CountReadsPerBin([test.bamFile1, test.bamFile2], 50, 4)
    >>> np.transpose(c.run())
    array([[ 0.,  0.,  1.,  1.],
           [ 0.,  1.,  1.,  2.]])
    """

    def __init__(self, bamFilesList, binLength=50, numberOfSamples=None, numberOfProcessors=1,
                 verbose=False, region=None,
                 bedFile=None, extendReads=False,
                 blackListFileName=None,
                 minMappingQuality=None,
                 ignoreDuplicates=False,
                 chrsToSkip=[],
                 stepSize=None,
                 center_read=False,
                 samFlag_include=None,
                 samFlag_exclude=None,
                 zerosToNans=False,
                 smoothLength=0,
                 minFragmentLength=0,
                 maxFragmentLength=0,
                 out_file_for_raw_data=None):

        self.bamFilesList = bamFilesList
        self.binLength = binLength
        self.numberOfSamples = numberOfSamples
        self.blackListFileName = blackListFileName

        if extendReads and len(bamFilesList):
            from deeptools.getFragmentAndReadSize import get_read_and_fragment_length
            frag_len_dict, read_len_dict = get_read_and_fragment_length(bamFilesList[0],
                                                                        return_lengths=False,
                                                                        blackListFileName=blackListFileName,
                                                                        numberOfProcessors=numberOfProcessors,
                                                                        verbose=verbose)
            if extendReads is True:
                # try to guess fragment length if the bam file contains paired end reads
                if frag_len_dict:
                    self.defaultFragmentLength = int(frag_len_dict['median'])
                else:
                    exit("*ERROR*: library is not paired-end. Please provide an extension length.")
                if verbose:
                    print(("Fragment length based on paired en data "
                          "estimated to be {}".format(frag_len_dict['median'])))

            elif extendReads < read_len_dict['median']:
                sys.stderr.write("*WARNING*: read extension is smaller than read length (read length = {}). "
                                 "Reads will not be extended.\n".format(int(read_len_dict['median'])))
                self.defaultFragmentLength = 'read length'

            elif extendReads > 2000:
                exit("*ERROR*: read extension must be smaller that 2000. Value give: {} ".format(extendReads))
            else:
                self.defaultFragmentLength = int(extendReads)

        else:
            self.defaultFragmentLength = 'read length'

        self.numberOfProcessors = numberOfProcessors
        self.verbose = verbose
        self.region = region
        self.bedFile = bedFile
        self.minMappingQuality = minMappingQuality
        self.ignoreDuplicates = ignoreDuplicates
        self.chrsToSkip = chrsToSkip
        self.stepSize = stepSize
        self.center_read = center_read
        self.samFlag_include = samFlag_include
        self.samFlag_exclude = samFlag_exclude
        self.minFragmentLength = minFragmentLength
        self.maxFragmentLength = maxFragmentLength
        self.zerosToNans = zerosToNans
        self.smoothLength = smoothLength

        if out_file_for_raw_data:
            self.save_data = True
            self.out_file_for_raw_data = out_file_for_raw_data
        else:
            self.save_data = False
            self.out_file_for_raw_data = None

        # check that wither numberOfSamples or stepSize are set
        if numberOfSamples is None and stepSize is None and bedFile is None:
            raise ValueError("either stepSize, numberOfSamples or bedFile have to be set")

        if self.defaultFragmentLength != 'read length':
            self.maxPairedFragmentLength = 4 * self.defaultFragmentLength
        else:
            self.maxPairedFragmentLength = 1000

    def run(self, allArgs=None):
        # Try to determine an optimal fraction of the genome (chunkSize) that is sent to
        # workers for analysis. If too short, too much time is spend loading the files
        # if too long, some processors end up free.
        # the following values are empirical

        bamFilesHandlers = [bamHandler.openBam(x) for x in self.bamFilesList]
        chromSizes, non_common = deeptools.utilities.getCommonChrNames(bamFilesHandlers, verbose=self.verbose)

        # skip chromosome in the list. This is usually for the
        # X chromosome which may have either one copy  in a male sample
        # or a mixture of male/female and is unreliable.
        # Also the skip may contain heterochromatic regions and
        # mitochondrial DNA
        if len(self.chrsToSkip):
            chromSizes = [x for x in chromSizes if x[0] not in self.chrsToSkip]

        chrNames, chrLengths = list(zip(*chromSizes))

        genomeSize = sum(chrLengths)
        if self.stepSize is None:
            if self.region is None:
                self.stepSize = max(int(float(genomeSize) / self.numberOfSamples), 1)
            else:
                # compute the step size, based on the number of samples
                # and the length of the region studied
                (chrom, start, end) = mapReduce.getUserRegion(chromSizes, self.region)[:3]
                self.stepSize = max(int(float(end - start) / self.numberOfSamples), 1)

        # number of samples is better if large
        if np.mean(chrLengths) < self.stepSize and self.bedFile is None:
            min_num_of_samples = int(genomeSize / np.mean(chrLengths))
            raise ValueError("numberOfSamples has to be bigger than {} ".format(min_num_of_samples))

        max_mapped = max([x.mapped for x in bamFilesHandlers])

        reads_per_bp = float(max_mapped) / genomeSize
        # chunkSize =  int(100 / ( reads_per_bp  * len(bamFilesList)) )

        chunkSize = int(self.stepSize * 1e3 / (reads_per_bp * len(bamFilesHandlers)))
        [bam_h.close() for bam_h in bamFilesHandlers]

        if self.verbose:
            print("step size is {}".format(self.stepSize))

        if self.region:
            # in case a region is used, append the tilesize
            self.region += ":{}".format(self.binLength)

        # Handle GTF options
        transcriptID, exonID, transcript_id_designator, keepExons = deeptools.utilities.gtfOptions(allArgs)

        # use map reduce to call countReadsInRegions_wrapper
        imap_res = mapReduce.mapReduce([],
                                       countReadsInRegions_wrapper,
                                       chromSizes,
                                       self_=self,
                                       genomeChunkLength=chunkSize,
                                       bedFile=self.bedFile,
                                       blackListFileName=self.blackListFileName,
                                       region=self.region,
                                       numberOfProcessors=self.numberOfProcessors,
                                       transcriptID=transcriptID,
                                       exonID=exonID,
                                       keepExons=keepExons,
                                       transcript_id_designator=transcript_id_designator)

        if self.out_file_for_raw_data:
            if len(non_common):
                sys.stderr.write("*Warning*\nThe resulting bed file does not contain information for "
                                 "the chromosomes that were not common between the bigwig files\n")

            # concatenate intermediary bedgraph files
            ofile = open(self.out_file_for_raw_data, "w")
            for _values, tempFileName in imap_res:
                if tempFileName:
                    # concatenate all intermediate tempfiles into one
                    _foo = open(tempFileName, 'r')
                    shutil.copyfileobj(_foo, ofile)
                    _foo.close()
                    os.remove(tempFileName)

            ofile.close()

        try:
            num_reads_per_bin = np.concatenate([x[0] for x in imap_res], axis=0)
            return num_reads_per_bin

        except ValueError:
            if self.bedFile:
                sys.exit('\nNo coverage values could be computed.\n\n'
                         'Please check that the chromosome names in the BED file are found on the bam files.\n\n'
                         'The valid chromosome names are:\n{}'.format(chrNames))
            else:
                sys.exit('\nNo coverage values could be computed.\n\nCheck that all bam files are valid and '
                         'contain mapped reads.')

    def count_reads_in_region(self, chrom, start, end, bed_regions_list=None):
        """Counts the reads in each bam file at each 'stepSize' position
        within the interval (start, end) for a window or bin of size binLength.

        The stepSize controls the distance between bins. For example,
        a step size of 20 and a bin size of 20 will create bins next to
        each other. If the step size is smaller than the bin size the
        bins will overlap.

        If a list of bedRegions is given, then the number of reads
        that overlaps with each region is counted.

        Parameters
        ----------
        chrom : str
            Chrom name
        start : int
            start coordinate
        end : int
            end coordinate
        bed_regions_list: list
            List of list of tuples of the form (start, end)
            corresponding to bed regions to be processed.
            If not bed file was passed to the object constructor
            then this list is empty.

        Returns
        -------
        numpy array
            The result is a numpy array that as rows each bin
            and as columns each bam file.


        Examples
        --------
        Initialize some useful values

        >>> test = Tester()
        >>> c = CountReadsPerBin([test.bamFile1, test.bamFile2], 25, 0, stepSize=50)

        The transpose is used to get better looking numbers. The first line
        corresponds to the number of reads per bin in the first bamfile.

        >>> _array, __ = c.count_reads_in_region(test.chrom, 0, 200)
        >>> _array
        array([[ 0.,  0.],
               [ 0.,  1.],
               [ 1.,  1.],
               [ 1.,  2.]])

        """

        if start > end:
            raise NameError("start %d bigger that end %d" % (start, end))

        if self.stepSize is None:
            raise ValueError("stepSize is not set!")
        # array to keep the read counts for the regions
        subnum_reads_per_bin = []

        start_time = time.time()

        bam_handlers = [bamHandler.openBam(bam) for bam in self.bamFilesList]

        blackList = None
        if self.blackListFileName is not None:
            blackList = GTF(self.blackListFileName)

        # A list of lists of tuples
        transcriptsToConsider = []
        if bed_regions_list is not None:
            transcriptsToConsider = [x[1] for x in bed_regions_list]
        else:
            if self.stepSize == self.binLength:
                transcriptsToConsider.append([(start, end, self.binLength)])
            else:
                for i in range(start, end, self.stepSize):
                    if i + self.binLength > end:
                        break
                    if blackList is not None and blackList.findOverlaps(chrom, i, i + self.binLength):
                        continue
                    transcriptsToConsider.append([(i, i + self.binLength)])

        if self.save_data:
            _file = open(deeptools.utilities.getTempFileName(suffix='.bed'), 'w+t')
            _file_name = _file.name
        else:
            _file_name = ''

        for bam in bam_handlers:
            for trans in transcriptsToConsider:
                tcov = self.get_coverage_of_region(bam, chrom, trans)
                if bed_regions_list is not None:
                    subnum_reads_per_bin.append(np.sum(tcov))
                else:
                    subnum_reads_per_bin.extend(tcov)

        subnum_reads_per_bin = np.concatenate([subnum_reads_per_bin]).reshape(-1, len(self.bamFilesList), order='F')

        if self.save_data:
            idx = 0
            for i, trans in enumerate(transcriptsToConsider):
                if len(trans[0]) != 3:
                    starts = ",".join([str(x[0]) for x in trans])
                    ends = ",".join([str(x[1]) for x in trans])
                    _file.write("\t".join([chrom, starts, ends]) + "\t")
                    _file.write("\t".join(["{}".format(x) for x in subnum_reads_per_bin[i, :]]) + "\n")
                else:
                    for exon in trans:
                        for startPos in range(exon[0], exon[1], exon[2]):
                            if idx >= subnum_reads_per_bin.shape[0]:
                                # At the end of chromosomes (or due to blacklisted regions), there are bins smaller than the bin size
                                # Counts there are added to the bin before them, but range() will still try to include them.
                                break
                            _file.write("{0}\t{1}\t{2}\t".format(chrom, startPos, startPos + exon[2]))
                            _file.write("\t".join(["{}".format(x) for x in subnum_reads_per_bin[idx, :]]) + "\n")
                            idx += 1
            _file.close()

        if self.verbose:
            endTime = time.time()
            rows = subnum_reads_per_bin.shape[0]
            print("%s countReadsInRegions_worker: processing %d "
                  "(%.1f per sec) @ %s:%s-%s" %
                  (multiprocessing.current_process().name,
                   rows, rows / (endTime - start_time), chrom, start, end))

        return subnum_reads_per_bin, _file_name

    def get_coverage_of_region(self, bamHandle, chrom, regions,
                               fragmentFromRead_func=None):
        """
        Returns a numpy array that corresponds to the number of reads
        that overlap with each tile.

        >>> test = Tester()
        >>> import pysam
        >>> c = CountReadsPerBin([], stepSize=1, extendReads=300)

        For this case the reads are length 36. The number of overlapping
        read fragments is 4 and 5 for the positions tested.

        >>> c.get_coverage_of_region(pysam.AlignmentFile(test.bamFile_PE), 'chr2',
        ... [(5000833, 5000834), (5000834, 5000835)])
        array([ 4.,  5.])

        In the following example a paired read is extended to the fragment length which is 100
        The first mate starts at 5000000 and the second at 5000064. Each mate is
        extended to the fragment length *independently*
        At position 500090-500100 one fragment  of length 100 overlap, and after position 5000101
        there should be zero reads.

        >>> c.zerosToNans = True
        >>> c.get_coverage_of_region(pysam.AlignmentFile(test.bamFile_PE), 'chr2',
        ... [(5000090, 5000100), (5000100, 5000110)])
        array([  1.,  nan])

        In the following  case the reads length is 50. Reads are not extended.

        >>> c.extendReads=False
        >>> c.get_coverage_of_region(pysam.AlignmentFile(test.bamFile2), '3R', [(148, 150), (150, 152), (152, 154)])
        array([ 1.,  2.,  2.])


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
            if chrom in bamHandle.references:
                reads = [r for r in bamHandle.fetch(chrom, regStart, regEnd)
                         if r.flag & 4 == 0]
            else:
                raise NameError("chromosome {} not found in bam file".format(chrom))

            prev_start_pos = None  # to store the start positions
            # of previous processed read pair
            for read in reads:
                if self.minMappingQuality and read.mapq < self.minMappingQuality:
                    continue

                # filter reads based on SAM flag
                if self.samFlag_include and read.flag & self.samFlag_include != self.samFlag_include:
                    continue
                if self.samFlag_exclude and read.flag & self.samFlag_exclude != 0:
                    continue

                # Fragment lengths
                if self.minFragmentLength > 0 and abs(read.template_length) < self.minFragmentLength:
                    continue
                if self.maxFragmentLength > 0 and abs(read.template_length) > self.maxFragmentLength:
                    continue

                # get rid of duplicate reads that have same position on each of the
                # pairs
                if self.ignoreDuplicates and prev_start_pos \
                        and prev_start_pos == (read.reference_start, read.pnext, read.is_reverse):
                    continue

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

                    sIdx = vector_start + max((fragmentStart - reg[0]) // tileSize, 0)
                    eIdx = vector_start + min(np.ceil(float(fragmentEnd - reg[0]) / tileSize).astype('int'), nRegBins)
                    if last_eIdx is not None:
                        sIdx = max(last_eIdx, sIdx)
                        if sIdx >= eIdx:
                            continue
                    coverages[sIdx:eIdx] += 1
                    last_eIdx = eIdx

                prev_start_pos = (read.reference_start, read.pnext, read.is_reverse)
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

    def getReadLength(self, read):
        return len(read)

    @staticmethod
    def is_proper_pair(read, maxPairedFragmentLength):
        """
        Checks if a read is proper pair meaning that both mates are facing each other and are in
        the same chromosome and are not to far away. The sam flag for proper pair can not
        always be trusted.
        :return: bool

        >>> import pysam
        >>> import os
        >>> from deeptools.countReadsPerBin import CountReadsPerBin as cr
        >>> root = os.path.dirname(os.path.abspath(__file__)) + "/test/test_data/"
        >>> bam = pysam.AlignmentFile("{}/test_proper_pair_filtering.bam".format(root))
        >>> iter = bam.fetch()
        >>> read = next(iter)
        >>> cr.is_proper_pair(read, 1000) # "keep" read
        True
        >>> cr.is_proper_pair(read, 200) # "keep" read, but maxPairedFragmentLength is too short
        False
        >>> read = next(iter)
        >>> cr.is_proper_pair(read, 1000) # "improper pair"
        False
        >>> read = next(iter)
        >>> cr.is_proper_pair(read, 1000) # "mismatch chr"
        False
        >>> read = next(iter)
        >>> cr.is_proper_pair(read, 1000) # "same orientation1"
        False
        >>> read = next(iter)
        >>> cr.is_proper_pair(read, 1000) # "same orientation2"
        False
        >>> read = next(iter)
        >>> cr.is_proper_pair(read, 1000) # "rev first"
        False
        >>> read = next(iter)
        >>> cr.is_proper_pair(read, 1000) # "rev first OK"
        True
        >>> read = next(iter)
        >>> cr.is_proper_pair(read, 1000) # "for first"
        False
        >>> read = next(iter)
        >>> cr.is_proper_pair(read, 1000) # "for first"
        True
        """
        if not read.is_proper_pair:
            return False
        if read.reference_id != read.next_reference_id:
            return False
        if abs(read.template_length) > maxPairedFragmentLength:
            return False
        # check that the mates face each other (inward)
        if read.is_reverse is read.mate_is_reverse:
            return False
        if read.is_reverse:
            if read.reference_start >= read.next_reference_start:
                return True
        else:
            if read.reference_start <= read.next_reference_start:
                return True
        return False

    def get_fragment_from_read(self, read):
        """Get read start and end position of a read.
        If given, the reads are extended as follows:
        If reads are paired end, each read mate is extended to match
        the fragment length, otherwise, a default fragment length
        is used. If reads are split (give by the CIGAR string) then
        the multiple positions of the read are returned.
        When reads are extended the cigar information is
        skipped.

        Parameters
        ----------
        read: pysam object.

        The following values are defined (for forward reads)::


                 |--          -- read.tlen --              --|
                 |-- read.alen --|
            -----|===============>------------<==============|----
                 |               |            |
            read.reference_start
                        read.reference_end  read.pnext

              and for reverse reads


                 |--             -- read.tlen --           --|
                                             |-- read.alen --|
            -----|===============>-----------<===============|----
                 |                           |               |
              read.pnext           read.reference_start  read.reference_end

        this is a sketch of a pair-end reads

        The function returns the fragment start and end, either
        using the paired end information (if available) or
        extending the read in the appropriate direction if this
        is single-end.

        Parameters
        ----------
        read : pysam read object


        Returns
        -------
        list of tuples
            [(fragment start, fragment end)]


        >>> test = Tester()
        >>> c = CountReadsPerBin([], 1, 1, 200, extendReads=True)
        >>> c.defaultFragmentLength=100
        >>> c.get_fragment_from_read(test.getRead("paired-forward"))
        [(5000000, 5000100)]
        >>> c.get_fragment_from_read(test.getRead("paired-reverse"))
        [(5000000, 5000100)]
        >>> c.defaultFragmentLength = 200
        >>> c.get_fragment_from_read(test.getRead("single-forward"))
        [(5001491, 5001691)]
        >>> c.get_fragment_from_read(test.getRead("single-reverse"))
        [(5001536, 5001736)]
        >>> c.defaultFragmentLength = 'read length'
        >>> c.get_fragment_from_read(test.getRead("single-forward"))
        [(5001491, 5001527)]
        >>> c.defaultFragmentLength = 'read length'
        >>> c.extendReads = False
        >>> c.get_fragment_from_read(test.getRead("paired-forward"))
        [(5000000, 5000036)]

        Tests for read centering.

        >>> c = CountReadsPerBin([], 1, 1, 200, extendReads=True, center_read=True)
        >>> c.defaultFragmentLength = 100
        >>> assert(c.get_fragment_from_read(test.getRead("paired-forward")) == [(5000032, 5000068)])
        >>> c.defaultFragmentLength = 200
        >>> assert(c.get_fragment_from_read(test.getRead("single-reverse")) == [(5001618, 5001654)])
        """
        # if no extension is needed, use pysam get_blocks
        # to identify start and end reference positions.
        # get_blocks return a list of start and end positions
        # based on the CIGAR if skipped regions are found.
        # E.g for a cigar of 40M260N22M
        # get blocks return two elements for the first 40 matches
        # and the for the last 22 matches.
        if self.defaultFragmentLength == 'read length':
            return read.get_blocks()

        else:
            if self.is_proper_pair(read, self.maxPairedFragmentLength):
                if read.is_reverse:
                    fragmentStart = read.next_reference_start
                    fragmentEnd = read.reference_end
                else:
                    fragmentStart = read.reference_start
                    # the end of the fragment is defined as
                    # the start of the forward read plus the insert length
                    fragmentEnd = read.reference_start + abs(read.template_length)

            # Extend using the default fragment length
            else:
                if read.is_reverse:
                    fragmentStart = read.reference_end - self.defaultFragmentLength
                    fragmentEnd = read.reference_end
                else:
                    fragmentStart = read.reference_start
                    fragmentEnd = read.reference_start + self.defaultFragmentLength

        if self.center_read:
            fragmentCenter = fragmentEnd - (fragmentEnd - fragmentStart) / 2
            fragmentStart = fragmentCenter - read.infer_query_length(always=False) / 2
            fragmentEnd = fragmentStart + read.infer_query_length(always=False)

        assert fragmentStart < fragmentEnd, "fragment start greater than fragment" \
                                            "end for read {}".format(read.query_name)
        return [(fragmentStart, fragmentEnd)]

    def getSmoothRange(self, tileIndex, tileSize, smoothRange, maxPosition):
        """
        Given a tile index position and a tile size (length), return the a new indices
        over a larger range, called the smoothRange.
        This region is centered in the tileIndex  an spans on both sizes
        to cover the smoothRange. The smoothRange is trimmed in case it is less
        than zero or greater than  maxPosition ::


             ---------------|==================|------------------
                        tileStart
                   |--------------------------------------|
                   |    <--      smoothRange     -->      |
                   |
             tileStart - (smoothRange-tileSize)/2

        Test for a smooth range that spans 3 tiles.

        Examples
        --------

        >>> c = CountReadsPerBin([], 1, 1, 1, 0)
        >>> c.getSmoothRange(5, 1, 3, 10)
        (4, 7)

        Test smooth range truncated on start.

        >>> c.getSmoothRange(0, 10, 30, 200)
        (0, 2)

        Test smooth range truncated on start.

        >>> c.getSmoothRange(1, 10, 30, 4)
        (0, 3)

        Test smooth range truncated on end.

        >>> c.getSmoothRange(5, 1, 3, 5)
        (4, 5)

        Test smooth range not multiple of tileSize.

        >>> c.getSmoothRange(5, 10, 24, 10)
        (4, 6)
        """
        smoothTiles = int(smoothRange / tileSize)
        if smoothTiles == 1:
            return (tileIndex, tileIndex + 1)

        smoothTilesSide = float(smoothTiles - 1) / 2
        smoothTilesLeft = int(np.ceil(smoothTilesSide))
        smoothTilesRight = int(np.floor(smoothTilesSide)) + 1

        indexStart = max(tileIndex - smoothTilesLeft, 0)
        indexEnd = min(maxPosition, tileIndex + smoothTilesRight)
        return (indexStart, indexEnd)


def remove_row_of_zeros(matrix):
    # remove rows containing all zeros or all nans
    _mat = np.nan_to_num(matrix)
    to_keep = _mat.sum(1) != 0
    return matrix[to_keep, :]


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
        self.root = os.path.dirname(os.path.abspath(__file__)) + "/test/test_data/"
        # self.root = "./test/test_data/"
        self.bamFile1 = self.root + "testA.bam"
        self.bamFile2 = self.root + "testB.bam"
        self.bamFile_PE = self.root + "test_paired2.bam"
        self.chrom = '3R'
        global debug
        debug = 0

    def getRead(self, readType):
        """ prepare arguments for test
        """
        bam = bamHandler.openBam(self.bamFile_PE)
        if readType == 'paired-reverse':
            read = [x for x in bam.fetch('chr2', 5000081, 5000082)][0]
        elif readType == 'single-forward':
            read = [x for x in bam.fetch('chr2', 5001491, 5001492)][0]
        elif readType == 'single-reverse':
            read = [x for x in bam.fetch('chr2', 5001700, 5001701)][0]
        else:  # by default a forward paired read is returned
            read = [x for x in bam.fetch('chr2', 5000027, 5000028)][0]
        return read
