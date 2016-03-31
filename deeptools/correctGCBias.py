#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import time
import subprocess
import sys

import twobitreader as twobit
import pysam
import multiprocessing
import numpy as np
import argparse

from scipy.stats import binom

from deeptools.utilities import getGC_content, tbitToBamChrName
from deeptools import writeBedGraph, parserCommon, mapReduce
from deeptools import utilities

old_settings = np.seterr(all='ignore')


def parse_arguments(args=None):
    parentParser = parserCommon.getParentArgParse(binSize=True, blackList=False)
    requiredArgs = getRequiredArgs()
    parser = argparse.ArgumentParser(
        parents=[requiredArgs, parentParser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='This tool corrects the GC-bias using the'
        ' method proposed by [Benjamini & Speed (2012). '
        'Nucleic Acids Research, 40(10)]. It will remove reads'
        ' from regions with too high coverage compared to the'
        ' expected values (typically GC-rich regions) and will'
        ' add reads to regions where too few reads are seen '
        '(typically AT-rich regions). '
        'The tool ``computeGCBias`` needs to be run first to generate the '
        'frequency table needed here.',
        usage='An example usage is:\n correctGCBias '
        '-b file.bam --effectiveGenomeSize 2150570000 -g mm9.2bit '
        '--GCbiasFrequenciesFile freq.txt -o gc_corrected.bam '
        '[options]',
        conflict_handler='resolve',
        add_help=False)
    return parser


def process_args(args=None):
    args = parse_arguments().parse_args(args)

    return args


def getRequiredArgs():
    parser = argparse.ArgumentParser(add_help=False)

    required = parser.add_argument_group('Required arguments')

    # define the arguments
    required.add_argument('--bamfile', '-b',
                          metavar='BAM file',
                          help='Sorted BAM file to correct.',
                          required=True)
    required.add_argument('--effectiveGenomeSize',
                          help='The effective genome size is the portion '
                          'of the genome that is mappable. Large fractions of '
                          'the genome are stretches of NNNN that should be '
                          'discarded. Also, if repetitive regions were not '
                          'included in the mapping of reads, the effective '
                          'genome size needs to be adjusted accordingly. '
                          'Common values are: mm9: 2150570000, '
                          'hg19:2451960000, dm3:121400000 and ce10:93260000. '
                          'See Table 2 of '
                          'http://www.plosone.org/article/info:doi/10.1371/journal.pone.0030377 '
                          'or http://www.nature.com/nbt/journal/v27/n1/fig_tab/nbt.1518_T1.html '
                          'for several effective genome sizes. This value is '
                          'needed to detect enriched regions that, if not '
                          'discarded, could bias the results.',
                          default=None,
                          type=int,
                          required=True)

    required.add_argument('--genome', '-g',
                          help='Genome in two bit format. Most genomes can be '
                          'found here: http://hgdownload.cse.ucsc.edu/gbdb/  '
                          'Search for the .2bit ending. Otherwise, fasta '
                          'files can be converted to 2bit using faToTwoBit '
                          'available here: '
                          'http://hgdownload.cse.ucsc.edu/admin/exe/',
                          metavar='two bit file',
                          required=True)

    required.add_argument('--GCbiasFrequenciesFile', '-freq',
                          help='Indicate the output file from '
                          'computeGCBias containing '
                          'the observed and expected read frequencies per GC-'
                          'content.',
                          type=argparse.FileType('r'),
                          metavar='FILE',
                          required=True)

    output = parser.add_argument_group('Output options')
    output.add_argument('--correctedFile', '-o',
                        help='Name of the corrected file. The ending will '
                        'be used to decide the output file format. The options '
                        'are ".bam", ".bw" for a bigWig file, ".bg" for a '
                        'bedGraph file.',
                        metavar='FILE',
                        type=argparse.FileType('w'),
                        required=True)

    # define the optional arguments
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")

    return parser


def getReadGCcontent(tbit, read, fragmentLength, chrNameBit):
    """
    The fragments for forward and reverse reads are defined as follows::

           |- read.pos       |- read.aend
        ---+=================>-----------------------+---------    Forward strand

           |-fragStart                               |-fragEnd

        ---+-----------------------<=================+---------    Reverse strand
                                   |-read.pos        |-read.aend

           |-----------------------------------------|
                            read.tlen

    """
    fragStart = None
    fragEnd = None

    if read.is_paired and read.is_proper_pair and abs(read.tlen) < 2 * fragmentLength:
        if read.is_reverse and read.tlen < 0:
            fragEnd = read.aend
            fragStart = read.aend + read.tlen
        elif read.tlen >= read.qlen:
            fragStart = read.pos
            fragEnd = read.pos + read.tlen

    if not fragStart:
        if read.is_reverse:
            fragEnd = read.aend
            fragStart = read.aend - fragmentLength
        else:
            fragStart = read.pos
            fragEnd = fragStart + fragmentLength
    fragStart = max(0, fragStart)
    try:
        gc = getGC_content(tbit[chrNameBit][fragStart:fragEnd], as_fraction=True)
    except Exception:
        return None
    if gc is None:
        return None

    # match the gc to the given fragmentLength
    gc = int(np.round(gc * fragmentLength))
    return gc


def writeCorrected_wrapper(args):
    return writeCorrected_worker(*args)


def writeCorrected_worker(chrNameBam, chrNameBit, start, end, step):
    r"""writes a bedgraph file containing the GC correction of
    a region from the genome

    >>> test = Tester()
    >>> tempFile = writeCorrected_worker(*test.testWriteCorrectedChunk())
    >>> open(tempFile, 'r').readlines()
    ['chr2L\t200\t225\t31.6\n', 'chr2L\t225\t250\t33.8\n', 'chr2L\t250\t275\t37.9\n', 'chr2L\t275\t300\t40.9\n']
    >>> os.remove(tempFile)
    """
    global R_gc
    fragmentLength = len(R_gc) - 1

    cvg_corr = np.zeros(end - start)

    i = 0

    tbit = twobit.TwoBitFile(global_vars['2bit'])
    bam = pysam.Samfile(global_vars['bam'])
    read_repetitions = 0
    removed_duplicated_reads = 0
    startTime = time.time()

    # caching seems to be faster
    # r.flag & 4 == 0 is to skip unmapped
    # reads that nevertheless are asigned
    # to a genomic position
    reads = [r for r in bam.fetch(chrNameBam, start, end)
             if r.flag & 4 == 0]

    bam.close()
    r_index = -1
    for read in reads:
        r_index += 1
        try:
            # calculate GC content of read fragment
            gc = getReadGCcontent(tbit, read, fragmentLength,
                                  chrNameBit)
        except Exception as detail:
            print(detail)
            """ this exception happens when the end of a
            chromosome is reached """
            continue
        if not gc:
            continue

        # is this read in the same orientation and position as the previous?
        if r_index > 0 and read.pos == reads[r_index - 1].pos and \
                read.is_reverse == reads[r_index - 1].is_reverse \
                and read.pnext == reads[r_index - 1].pnext:
            read_repetitions += 1
            if read_repetitions >= global_vars['max_dup_gc'][gc]:
                removed_duplicated_reads += 1
                continue
        else:
            read_repetitions = 0

        try:
            fragmentStart, fragmentEnd = getFragmentFromRead(read, fragmentLength, extendPairedEnds=True)
            vectorStart = max(fragmentStart - start, 0)
            vectorEnd = min(fragmentEnd - start, end - start)
        except TypeError:
            # the get_fragment_from_read functions returns None in some cases.
            # Those cases are to be skipped, hence the continue line.
            continue

        cvg_corr[vectorStart:vectorEnd] += float(1) / R_gc[gc]
        i += 1
    if debug:
        endTime = time.time()
        print("{}, processing {} ({:.1f} per sec) ")
        "reads @ {}:{}-{}".format(multiprocessing.current_process().name,
                                  i, i / (endTime - startTime),
                                  chrNameBit, start, end)

    if i == 0:
        return None

    _file = open(utilities.getTempFileName(suffix='.bg'), 'w')
    # save in bedgraph format
    for bin in range(0, len(cvg_corr), step):
        value = np.mean(cvg_corr[bin:min(bin + step, end)])
        if value > 0:
            writeStart = start + bin
            writeEnd = min(start + bin + step, end)
            _file.write("%s\t%d\t%d\t%.1f\n" % (chrNameBit, writeStart,
                                                writeEnd, value))

    tempFileName = _file.name
    _file.close()
    return tempFileName


def numCopiesOfRead(value):
    """
    Based int he R_gc value, decides
    whether to keep, duplicate, triplicate or delete the read.
    It returns an integer, that tells the number of copies of the read
    that should be keep.
    >>> np.random.seed(1)
    >>> numCopiesOfRead(0.8)
    1
    >>> numCopiesOfRead(2.5)
    2
    >>> numCopiesOfRead(None)
    1
    """
    copies = 1
    if value:
        copies = int(value) + (1 if np.random.rand() < value % 1 else 0)
    return copies


def writeCorrectedSam_wrapper(args):
    return writeCorrectedSam_worker(*args)


def writeCorrectedSam_worker(chrNameBam, chrNameBit, start, end,
                             step=None,
                             tag_but_not_change_number=False,
                             verbose=True):
    r"""
    Writes a BAM file, deleting and adding some reads in order to compensate
    for the GC bias. **This is a stochastic method.**
    >>> np.random.seed(1)
    >>> test = Tester()
    >>> args = test.testWriteCorrectedSam()
    >>> tempFile = writeCorrectedSam_worker(*args, \
    ... tag_but_not_change_number=True, verbose=False)
    >>> try:
    ...     import StringIO
    ... except ImportError:
    ...     from io import StringIO
    >>> ostdout = sys.stdout
    >>> import tempfile
    >>> sys.stdout = tempfile.TemporaryFile()
    >>> idx = pysam.index(tempFile)
    >>> sys.stdout = ostdout
    >>> bam = pysam.Samfile(tempFile)
    >>> [dict(r.tags)['YN'] for r in bam.fetch(args[0], 200, 250)]
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1]
    >>> res = os.remove(tempFile)
    >>> res = os.remove(tempFile+".bai")
    >>> tempFile = \
    ... writeCorrectedSam_worker(*test.testWriteCorrectedSam_paired(),\
    ... tag_but_not_change_number=True, verbose=False)
    >>> sys.stdout = tempfile.TemporaryFile()
    >>> idx = pysam.index(tempFile)
    >>> sys.stdout = ostdout
    >>> bam = pysam.Samfile(tempFile)
    >>> [dict(r.tags)['YN'] for r in bam.fetch('chr2L', 0, 50)]
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    >>> res = os.remove(tempFile)
    >>> res = os.remove(tempFile+".bai")
    """
    global R_gc
    fragmentLength = len(R_gc) - 1

    if verbose:
        print("Sam for %s %s %s " % (chrNameBit, start, end))
    i = 0

    tbit = twobit.TwoBitFile(global_vars['2bit'])

    bam = pysam.Samfile(global_vars['bam'])
    tempFileName = utilities.getTempFileName(suffix='.bam')

    outfile = pysam.Samfile(tempFileName, 'wb', template=bam)
    startTime = time.time()
    matePairs = {}
    read_repetitions = 0
    removed_duplicated_reads = 0
    # cache data
    # r.flag & 4 == 0 is to filter unmapped reads that
    # have a genomic position
    reads = [r for r in bam.fetch(chrNameBam, start, end)
             if r.pos > start and r.flag & 4 == 0]

    r_index = -1
    for read in reads:
        r_index += 1
        copies = None
        gc = None

        # check if a mate has already been procesed
        # to apply the same correction
        try:
            copies = matePairs[read.qname]['copies']
            gc = matePairs[read.qname]['gc']
            del(matePairs[read.qname])
        except:
            # this exception happens when a mate is
            # not present. This could
            # happen because of removal of the mate
            # by some filtering
            gc = getReadGCcontent(tbit, read, fragmentLength,
                                  chrNameBit)
            if gc:
                copies = numCopiesOfRead(float(1) / R_gc[gc])
            else:
                copies = 1
        # is this read in the same orientation and position as the previous?
        if gc and r_index > 0 and read.pos == reads[r_index - 1].pos \
                and read.is_reverse == reads[r_index - 1].is_reverse \
                and read.pnext == reads[r_index - 1].pnext:
            read_repetitions += 1
            if read_repetitions >= global_vars['max_dup_gc'][gc]:
                copies = 0  # in other words do not take into account this read
                removed_duplicated_reads += 1
        else:
            read_repetitions = 0

        readName = read.qname
        # Each tag is a tuple of (tag name, value, type)
        # Note that get_tags() returns ord(type) rather than type and this must
        # be fixed!
        # It turns out that the "with_value_type" option only started working in
        # pysam-0.8.4, so we can't reliably add tags on earlier versions without
        # potentially creating BAM files that break HTSJDK/IGV/etc.

        readTag = read.get_tags(with_value_type=True)
        replace_tags = False
        if len(readTag) > 0:
            if len(readTag[0]) == 3:
                if type(readTag[2]) is int:
                    readTag = [(x[0], x[1], chr(x[2])) for x in readTag]
                replace_tags = True
        else:
            replace_tags = True

        if gc:
            GC = int(100 * np.round(float(gc) / fragmentLength,
                                    decimals=2))
            readTag.append(
                ('YC', float(round(float(1) / R_gc[gc], 2)), "f"))
            readTag.append(('YN', copies, "i"))
        else:
            GC = -1

        readTag.append(('YG', GC, "i"))
        if replace_tags:
            read.set_tags(readTag)

        if read.is_paired and read.is_proper_pair \
                and not read.mate_is_unmapped \
                and not read.is_reverse:
            matePairs[readName] = {'copies': copies,
                                   'gc': gc}

        """
        outfile.write(read)
        """
        if tag_but_not_change_number:
            outfile.write(read)
            continue

        for numCop in range(1, copies + 1):
            # the read has to be renamed such that newly
            # formed pairs will match
            if numCop > 1:
                read.qname = readName + "_%d" % (numCop)
            outfile.write(read)

        if verbose:
            if i % 500000 == 0 and i > 0:
                endTime = time.time()
                print("{},  processing {} ({:.1f} per sec) reads "
                      "@ {}:{}-{}".format(multiprocessing.current_process().name,
                                          i, i / (endTime - startTime),
                                          chrNameBit, start, end))
        i += 1

    outfile.close()
    if verbose:
        endTime = time.time()
        print("{},  processing {} ({:.1f} per sec) reads "
              "@ {}:{}-{}".format(multiprocessing.current_process().name,
                                  i, i / (endTime - startTime),
                                  chrNameBit, start, end))
        percentage = float(removed_duplicated_reads) * 100 / len(reads) \
            if len(reads) > 0 else 0
        print("duplicated reads removed %d of %d (%.2f) " %
              (removed_duplicated_reads, len(reads), percentage))

    return tempFileName


def getFragmentFromRead(read, defaultFragmentLength, extendPairedEnds=True):
    """
    The read has to be pysam object.

    The following values are defined (for forward reads)::


             |--          -- read.tlen --              --|
             |-- read.alen --|
        -----|===============>------------<==============|----
             |               |            |
          read.pos      read.aend      read.pnext


          and for reverse reads


             |--             -- read.tlen --           --|
                                         |-- read.alen --|
        -----|===============>-----------<===============|----
             |                           |               |
          read.pnext                   read.pos      read.aend

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
    tuple
        (fragment start, fragment end)

    """
    # convert reads to fragments

    # this option indicates that the paired ends correspond
    # to the fragment ends
    # condition read.tlen < maxPairedFragmentLength is added to avoid read pairs
    # that span thousands of base pairs

    if extendPairedEnds is True and read.is_paired and 0 < abs(read.tlen) < 1000:
        if read.is_reverse:
            fragmentStart = read.pnext
            fragmentEnd = read.aend
        else:
            fragmentStart = read.pos
            # the end of the fragment is defined as
            # the start of the forward read plus the insert length
            fragmentEnd = read.pos + read.tlen
    else:
        if defaultFragmentLength <= read.aend - read.pos:
            fragmentStart = read.pos
            fragmentEnd = read.aend
        else:
            if read.is_reverse:
                fragmentStart = read.aend - defaultFragmentLength
                fragmentEnd = read.aend
            else:
                fragmentStart = read.pos
                fragmentEnd = read.pos + defaultFragmentLength

    return fragmentStart, fragmentEnd


def run_shell_command(command):
    """
    Runs the given shell command. Report
    any errors found.
    """
    try:
        subprocess.check_call(command, shell=True)

    except subprocess.CalledProcessError as error:
        sys.stderr.write('Error{}\n'.format(error))
        exit(1)
    except Exception as error:
        sys.stderr.write('Error: {}\n'.format(error))
        exit(1)


def main(args=None):
    args = process_args(args)
    global F_gc, N_gc, R_gc

    data = np.loadtxt(args.GCbiasFrequenciesFile.name)

    F_gc = data[:, 0]
    N_gc = data[:, 1]
    R_gc = data[:, 2]

    global global_vars
    global_vars = {}
    global_vars['2bit'] = args.genome
    global_vars['bam'] = args.bamfile

    # compute the probability to find more than one read (a redundant read)
    # at a certain position based on the gc of the read fragment
    # the binomial function is used for that
    max_dup_gc = [binom.isf(1e-7, F_gc[x], 1.0 / N_gc[x])
                  if F_gc[x] > 0 and N_gc[x] > 0 else 1
                  for x in range(len(F_gc))]

    global_vars['max_dup_gc'] = max_dup_gc

    tbit = twobit.TwoBitFile(global_vars['2bit'])
    bam = pysam.Samfile(global_vars['bam'])

    global_vars['genome_size'] = sum(tbit.sequence_sizes().values())
    global_vars['total_reads'] = bam.mapped
    global_vars['reads_per_bp'] = \
        float(global_vars['total_reads']) / args.effectiveGenomeSize

    # apply correction
    print("applying correction")
    # divide the genome in fragments containing about 4e5 reads.
    # This amount of reads takes about 20 seconds
    # to process per core (48 cores, 256 Gb memory)
    chunkSize = int(4e5 / global_vars['reads_per_bp'])

    # chromSizes: list of tuples
    chromSizes = [(bam.references[i], bam.lengths[i])
                  for i in range(len(bam.references))]

    regionStart = 0
    if args.region:
        chromSizes, regionStart, regionEnd, chunkSize = \
            mapReduce.getUserRegion(chromSizes, args.region,
                                    max_chunk_size=chunkSize)

    print("genome partition size for multiprocessing: {}".format(chunkSize))
    print("using region {}".format(args.region))
    mp_args = []
    bedGraphStep = args.binSize
    chrNameBitToBam = tbitToBamChrName(list(tbit.sequence_sizes().keys()), bam.references)
    chrNameBamToBit = dict([(v, k) for k, v in chrNameBitToBam.items()])
    print(chrNameBitToBam, chrNameBamToBit)
    c = 1
    for chrom, size in chromSizes:
        start = 0 if regionStart == 0 else regionStart
        for i in range(start, size, chunkSize):
            try:
                chrNameBamToBit[chrom]
            except KeyError:
                print("no sequence information for ")
                "chromosome {} in 2bit file".format(chrom)
                print("Reads in this chromosome will be skipped")
                continue
            length = min(size, i + chunkSize)
            mp_args.append((chrom, chrNameBamToBit[chrom], i, length,
                            bedGraphStep))
            c += 1

    pool = multiprocessing.Pool(args.numberOfProcessors)

    if args.correctedFile.name.endswith('bam'):
        if len(mp_args) > 1 and args.numberOfProcessors > 1:
            print(("using {} processors for {} "
                   "number of tasks".format(args.numberOfProcessors,
                                            len(mp_args))))

            res = pool.map_async(
                writeCorrectedSam_wrapper, mp_args).get(9999999)
        else:
            res = list(map(writeCorrectedSam_wrapper, mp_args))

        if len(res) == 1:
            command = "cp {} {}".format(res[0], args.correctedFile.name)
            run_shell_command(command)
        else:
            print("concatenating (sorted) intermediate BAMs")
            header = pysam.Samfile(res[0])
            of = pysam.Samfile(args.correctedFile.name, "wb", template=header)
            header.close()
            for f in res:
                f = pysam.Samfile(f)
                for e in f.fetch(until_eof=True):
                    of.write(e)
                f.close()
            of.close()

        print("indexing BAM")
        pysam.index(args.correctedFile.name)

        for tempFileName in res:
            os.remove(tempFileName)

    if args.correctedFile.name.endswith('bg') or \
            args.correctedFile.name.endswith('bw'):

        _temp_bg_file_name = utilities.getTempFileName(suffix='_all.bg')
        if len(mp_args) > 1 and args.numberOfProcessors > 1:

            res = pool.map_async(writeCorrected_wrapper, mp_args).get(9999999)
        else:
            res = list(map(writeCorrected_wrapper, mp_args))

        # concatenate intermediary bedgraph files
        _temp_bg_file = open(_temp_bg_file_name, 'w')
        for tempFileName in res:
            if tempFileName:
                # concatenate all intermediate tempfiles into one
                # bedgraph file
                shutil.copyfileobj(open(tempFileName, 'rb'), _temp_bg_file)
                os.remove(tempFileName)
        _temp_bg_file.close()
        args.correctedFile.close()

        if args.correctedFile.name.endswith('bg'):
            shutil.move(_temp_bg_file_name, args.correctedFile.name)

        else:
            chromSizes = [(k, v) for k, v in tbit.sequence_sizes().items()]
            writeBedGraph.bedGraphToBigWig(chromSizes, _temp_bg_file_name,
                                           args.correctedFile.name)
            os.remove(_temp_bg_file)


class Tester():
    def __init__(self):
        import os
        self.root = os.path.dirname(os.path.abspath(__file__)) + "/test/test_corrGC/"
        self.tbitFile = self.root + "sequence.2bit"
        self.bamFile = self.root + "test.bam"
        self.chrNameBam = '2L'
        self.chrNameBit = 'chr2L'
        bam = pysam.Samfile(self.bamFile)
        tbit = twobit.TwoBitFile(self.tbitFile)
        global debug
        debug = 0
        global global_vars
        global_vars = {'2bit': self.tbitFile,
                       'bam': self.bamFile,
                       'filter_out': None,
                       'extra_sampling_file': None,
                       'max_reads': 5,
                       'min_reads': 0,
                       'min_reads': 0,
                       'reads_per_bp': 0.3,
                       'total_reads': bam.mapped,
                       'genome_size': sum(tbit.sequence_sizes().values())}

    def testWriteCorrectedChunk(self):
        """ prepare arguments for test
        """
        global R_gc, R_gc_min, R_gc_max
        R_gc = np.loadtxt(self.root + "R_gc_paired.txt")

        global_vars['max_dup_gc'] = np.ones(301)

        start = 200
        end = 300
        bedGraphStep = 25
        return (self.chrNameBam,
                self.chrNameBit, start, end, bedGraphStep)

    def testWriteCorrectedSam(self):
        """ prepare arguments for test
        """
        global R_gc, R_gc_min, R_gc_max
        R_gc = np.loadtxt(self.root + "R_gc_paired.txt")

        global_vars['max_dup_gc'] = np.ones(301)

        start = 200
        end = 250
        return (self.chrNameBam,
                self.chrNameBit, start, end)

    def testWriteCorrectedSam_paired(self):
        """ prepare arguments for test.
        """
        global R_gc, R_gc_min, R_gc_max
        R_gc = np.loadtxt(self.root + "R_gc_paired.txt")

        start = 0
        end = 500
        global global_vars
        global_vars['bam'] = self.root + "paired.bam"
        return 'chr2L', 'chr2L', start, end


if __name__ == "__main__":
    main()
