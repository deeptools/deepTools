import pyBigWig
import numpy as np
import os
import sys
import shutil
import warnings

# deepTools packages
import deeptools.mapReduce as mapReduce
import deeptools.utilities
# debug = 0

old_settings = np.seterr(all='ignore')


def countReadsInRegions_wrapper(args):
    # Using arguments unpacking!
    return countFragmentsInRegions_worker(*args)


def countFragmentsInRegions_worker(chrom, start, end,
                                   bigWigFiles,
                                   stepSize, binLength,
                                   save_data,
                                   bedRegions=None
                                   ):
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
    >>> np.transpose(countFragmentsInRegions_worker(test.chrom, 0, 200, [test.bwFile1, test.bwFile2], 50, 25, False)[0])
    array([[ 1.,  1.,  2.,  2.],
           [ 1.,  1.,  1.,  3.]])

    >>> np.transpose(countFragmentsInRegions_worker(test.chrom, 0, 200, [test.bwFile1, test.bwFile2], 200, 200, False)[0])
    array([[ 1.5],
           [ 1.5]])

    BED regions:
    >>> bedRegions = [[test.chrom, [(45, 55)]], [test.chrom, [(95, 105)]], [test.chrom, [(145, 155)]]]
    >>> np.transpose(countFragmentsInRegions_worker(test.chrom, 0, 200,[test.bwFile1, test.bwFile2], 200, 200, False,
    ... bedRegions=bedRegions)[0])
    array([[ 1. ,  1.5,  2. ],
           [ 1. ,  1. ,  2. ]])
    """
    assert start < end, "start {} bigger that end {}".format(start, end)

    # array to keep the scores for the regions
    sub_score_per_bin = []

    rows = 0

    bigwig_handlers = [pyBigWig.open(bw) for bw in bigWigFiles]

    regions_to_consider = []
    if bedRegions:
        for reg in bedRegions:
            regs = []
            for exon in reg[1]:
                regs.append((exon[0], exon[1]))
            regions_to_consider.append(regs)
    else:
        for i in range(start, end, stepSize):
            if (i + binLength) > end:
                regions_to_consider.append([(i, end)])  # last bin (may be smaller)
            else:
                regions_to_consider.append([(i, i + binLength)])

    if save_data:
        _file = open(deeptools.utilities.getTempFileName(suffix='.bed'), 'w+t')
        _file_name = _file.name
    else:
        _file_name = ''
    warnings.simplefilter("default")
    i = 0
    for reg in regions_to_consider:
        avgReadsArray = []
        i += 1

        for idx, bwh in enumerate(bigwig_handlers):
            if chrom not in list(bwh.chroms().keys()):
                unmod_name = chrom
                if chrom.startswith('chr'):
                    # remove the chr part from chromosome name
                    chrom = chrom[3:]
                else:
                    # prefix with 'chr' the chromosome name
                    chrom = 'chr' + chrom
                if chrom not in list(bwh.chroms().keys()):
                    exit('Chromosome name {} not found in bigwig file\n {}\n'.format(unmod_name, bigWigFiles[idx]))

            weights = []
            scores = []
            for exon in reg:
                weights.append(exon[1] - exon[0])
                score = bwh.stats(chrom, exon[0], exon[1])

                if score is None or score == [None] or np.isnan(score[0]):
                    score = [np.nan]
                scores.extend(score)
            avgReadsArray.append(np.average(scores, weights=weights))  # mean of fragment coverage for region

        sub_score_per_bin.extend(avgReadsArray)
        rows += 1
        if save_data:
            starts = []
            ends = []
            for exon in reg:
                starts.append(str(exon[0]))
                ends.append(str(exon[1]))
            starts = ",".join(starts)
            ends = ",".join(ends)
            _file.write("\t".join(map(str, [chrom, starts, ends])) + "\t")
            _file.write("\t".join(["{}".format(x) for x in avgReadsArray]) + "\n")

    if save_data:
        _file.close()
    warnings.resetwarnings()

    # the output is a matrix having as many rows as the variable 'row'
    # and as many columns as bigwig files. The rows correspond to
    # each of the regions processed by the worker.
    # np.array([[score1_1, score1_2],
    #           [score2_1, score2_2]]
    return np.array(sub_score_per_bin).reshape(rows, len(bigWigFiles)), _file_name


def getChromSizes(bigwigFilesList):
    """
    Get chromosome sizes from bigWig file with pyBigWig

    Test dataset with two samples covering 200 bp.
    >>> test = Tester()

    Chromosome name(s) and size(s).
    >>> assert(getChromSizes([test.bwFile1, test.bwFile2]) == ([('3R', 200)], set([])))
    """
    # check that the path to USCS bedGraphToBigWig as set in the config
    # is installed and is executable.

    def print_chr_names_and_size(chr_set):
        sys.stderr.write("chromosome\tlength\n")
        for name, size in chr_set:
            sys.stderr.write("{0:>15}\t{1:>10}\n".format(name, size))

    bigwigFilesList = bigwigFilesList[:]

    common_chr = set(pyBigWig.open(bigwigFilesList.pop()).chroms().items())
    non_common_chr = set()
    for bw in bigwigFilesList:
        _names_and_size = set(pyBigWig.open(bw).chroms().items())
        if len(common_chr & _names_and_size) == 0:
            #  try to add remove 'chr' from the chromosme name
            _corr_names_size = set()
            for chrom_name, size in _names_and_size:
                if chrom_name.startswith('chr'):
                    _corr_names_size.add((chrom_name[3:], size))
                else:
                    _corr_names_size.add(('chr' + chrom_name, size))
            if len(common_chr & _corr_names_size) == 0:
                message = "No common chromosomes found. Are the bigwig files " \
                          "from the same species and same assemblies?\n"
                sys.stderr.write(message)
                print_chr_names_and_size(common_chr)

                sys.stderr.write("\nand the following is the list of the unmatched chromosome and chromosome\n"
                                 "lengths from file\n{}\n".format(bw))
                print_chr_names_and_size(_names_and_size)
                exit(1)
            else:
                _names_and_size = _corr_names_size

        non_common_chr |= common_chr ^ _names_and_size
        common_chr = common_chr & _names_and_size

    if len(non_common_chr) > 0:
        sys.stderr.write("\nThe following chromosome names did not match between the the bigwig files\n")
        print_chr_names_and_size(non_common_chr)

    # get the list of common chromosome names and sizes
    return sorted(common_chr), non_common_chr


def getScorePerBin(bigWigFiles, binLength,
                   numberOfProcessors=1,
                   verbose=False, region=None,
                   bedFile=None,
                   blackListFileName=None,
                   stepSize=None,
                   chrsToSkip=[],
                   out_file_for_raw_data=None,
                   allArgs=None):
    """
    This function returns a matrix containing scores (median) for the coverage
    of fragments within a region. Each row corresponds to a sampled region.
    Likewise, each column corresponds to a bigwig file.

    Test dataset with two samples covering 200 bp.
    >>> test = Tester()
    >>> np.transpose(getScorePerBin([test.bwFile1, test.bwFile2], 50, 3))
    array([[ 1.,  1.,  2.,  2.],
           [ 1.,  1.,  1.,  3.]])

    """

    # Try to determine an optimal fraction of the genome (chunkSize)
    # that is sent to workers for analysis. If too short, too much time
    # is spent loading the files
    # if too long, some processors end up free.
    # the following is a heuristic

    # get list of common chromosome names and sizes
    chrom_sizes, non_common = getChromSizes(bigWigFiles)
    # skip chromosome in the list. This is usually for the
    # X chromosome which may have either one copy  in a male sample
    # or a mixture of male/female and is unreliable.
    # Also the skip may contain heterochromatic regions and
    # mitochondrial DNA
    if chrsToSkip and len(chrsToSkip):
        chrom_sizes = [x for x in chrom_sizes if x[0] not in chrsToSkip]

    chrnames, chrlengths = list(zip(*chrom_sizes))
    if stepSize is None:
        stepSize = binLength  # for adjacent bins

    # set chunksize based on number of processors used
    chunkSize = max(sum(chrlengths) / numberOfProcessors, int(1e6))
    # make chunkSize multiple of binLength
    chunkSize -= chunkSize % binLength
    if verbose:
        print("step size is {}".format(stepSize))

    if region:
        # in case a region is used, append the tilesize
        region += ":{}".format(binLength)
    # mapReduce( (staticArgs), func, chromSize, etc. )
    if out_file_for_raw_data:
        save_file = True
    else:
        save_file = False

    # Handle GTF options
    transcriptID, exonID, transcript_id_designator, keepExons = deeptools.utilities.gtfOptions(allArgs)

    imap_res = mapReduce.mapReduce((bigWigFiles, stepSize, binLength, save_file),
                                   countReadsInRegions_wrapper,
                                   chrom_sizes,
                                   genomeChunkLength=chunkSize,
                                   bedFile=bedFile,
                                   blackListFileName=blackListFileName,
                                   region=region,
                                   numberOfProcessors=numberOfProcessors,
                                   transcriptID=transcriptID,
                                   exonID=exonID,
                                   keepExons=keepExons,
                                   transcript_id_designator=transcript_id_designator)

    if out_file_for_raw_data:
        if len(non_common):
            sys.stderr.write("*Warning*\nThe resulting bed file does not contain information for "
                             "the chromosomes that were not common between the bigwig files\n")

        # concatenate intermediary bedgraph files
        ofile = open(out_file_for_raw_data, "w")
        for _values, tempFileName in imap_res:
            if tempFileName:
                # concatenate all intermediate tempfiles into one
                f = open(tempFileName, 'r')
                shutil.copyfileobj(f, ofile)
                f.close()
                os.remove(tempFileName)

        ofile.close()

    # the matrix scores are in the first element of each of the entries in imap_res
    score_per_bin = np.concatenate([x[0] for x in imap_res], axis=0)
    return score_per_bin


class Tester(object):

    def __init__(self):
        """
        The the two bigWig files are as follows:
        $ cat /tmp/testA.bg
        3R      0       100     1
        3R      100     200     2

        $ cat /tmp/testB.bg
        3R      0       150     1
        3R      150     200     3

        They cover 200 bp:

              0              50              100            150            200
              |------------------------------------------------------------|
            A  111111111111111111111111111111122222222222222222222222222222


            B  111111111111111111111111111111111111111111111333333333333333

        """

        self.root = os.path.dirname(os.path.abspath(__file__)) + "/test/test_data/"
        self.bwFile1 = self.root + "testA.bw"
        self.bwFile2 = self.root + "testB.bw"
        self.bwFile_PE = self.root + "test_paired2.bw"
        self.chrom = '3R'
        # global debug
        # debug = 0
