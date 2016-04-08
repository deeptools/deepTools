#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time

import multiprocessing
import numpy as np
import argparse
from scipy.stats import poisson
import twobitreader as twobit

from deeptoolsintervals import GTF
from deeptools.utilities import getGC_content, tbitToBamChrName
from deeptools import parserCommon, mapReduce
from deeptools.getFragmentAndReadSize import get_read_and_fragment_length
from deeptools import bamHandler

debug = 0
old_settings = np.seterr(all='ignore')


def parse_arguments(args=None):
    parentParser = parserCommon.getParentArgParse(binSize=False, blackList=True)
    requiredArgs = getRequiredArgs()
    parser = argparse.ArgumentParser(
        parents=[requiredArgs, parentParser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Computes the GC-bias using Benjamini\'s method '
        '[Benjamini & Speed (2012). Nucleic Acids Research, 40(10). doi: 10.1093/nar/gks001]. '
        'The GC-bias is visualized and the resulting table can be used to'
        'correct the bias with `correctGCBias`.',
        usage='\n computeGCBias '
        '-b file.bam --effectiveGenomeSize 2150570000 -g mm9.2bit -l 200 --GCbiasFrequenciesFile freq.txt [options]',
        conflict_handler='resolve',
        add_help=False)

    return parser


def getRequiredArgs():
    parser = argparse.ArgumentParser(add_help=False)

    required = parser.add_argument_group('Required arguments')

    required.add_argument('--bamfile', '-b',
                          metavar='bam file',
                          help='Sorted BAM file. ',
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
                          'discarded can bias the results.',
                          default=None,
                          type=int,
                          required=True)

    required.add_argument('--genome', '-g',
                          help='Genome in two bit format. Most genomes can be '
                          'found here: http://hgdownload.cse.ucsc.edu/gbdb/ '
                          'Search for the .2bit ending. Otherwise, fasta '
                          'files can be converted to 2bit using the UCSC '
                          'programm called faToTwoBit available for different '
                          'plattforms at '
                          'http://hgdownload.cse.ucsc.edu/admin/exe/',
                          metavar='2bit FILE',
                          required=True)

    required.add_argument('--fragmentLength', '-l',
                          help='Fragment length used for the sequencing. If '
                          'paired-end reads are used, the fragment length is '
                          'computed based from the bam file',
                          type=int,
                          required=True)

    # define the optional arguments
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")

    optional.add_argument('--sampleSize',
                          default=5e7,
                          help='Number of sampling points to be considered.',
                          type=int)

    optional.add_argument('--extraSampling',
                          help='BED file containing genomic regions for which '
                          'extra sampling is required because they are '
                          'underrepresented in the genome.',
                          type=argparse.FileType('r'),
                          metavar='BED file')

    group = parser.add_argument_group('Output options')

    group.add_argument('--GCbiasFrequenciesFile', '-freq',
                       help='Path to save the file containing '
                       'the observed and expected read frequencies per %%GC-'
                       'content. This file is needed to run the '
                       'correctGCBias tool. This is a text file.',
                       type=argparse.FileType('w'),
                       metavar='FILE',
                       required=True)

    plot = parser.add_argument_group('Diagnostic plot options')

    plot.add_argument('--biasPlot',
                      metavar='FILE NAME',
                      help='If given, a diagnostic image summarizing '
                      'the GC-bias will be saved.')

    plot.add_argument('--regionSize',
                      metavar='INT',
                      type=int,
                      default=300,
                      help='To plot the reads per %%GC over a region'
                      'the size of the region is required. By default, '
                      'the bin size is set to 300 bases, which is close to the '
                      'standard fragment size for Illumina machines. However, '
                      'if the depth of sequencing is low, a larger bin size '
                      'will be required, otherwise many bins will not '
                      'overlap with any read')

    group.add_argument('--plotFileFormat',
                       metavar='',
                       help='image format type. If given, this '
                       'option overrides the '
                       'image format based on the plotFile ending. '
                       'The available options are: "png", '
                       '"eps", "pdf" and "svg"',
                       choices=['png', 'pdf', 'svg', 'eps'])

    return parser


def getPositionsToSample(chrom, start, end, stepSize):
    """
    check if the region submitted to the worker
    overlaps with the region to take extra effort to sample.
    If that is the case, the regions to sample array is
    increased to match each of the positions in the extra
    effort region sampled at the same stepSize along the interval.

    If a filter out tree is given, then from positions to sample
    those regions are cleaned
    """
    positions_to_sample = np.arange(start, end, stepSize)

    if global_vars['filter_out']:
        filter_out_tree = GTF(global_vars['filter_out'])
    else:
        filter_out_tree = None

    if global_vars['extra_sampling_file']:
        extra_tree = GTF(global_vars['extra_sampling_file'])
    else:
        extra_tree = None

    if extra_tree:
        orig_len = len(positions_to_sample)
        try:
            extra_match = extra_tree.findOverlaps(chrom, start, end)
        except KeyError:
            extra_match = []

        if len(extra_match) > 0:
            for intval in extra_match:
                positions_to_sample = np.append(positions_to_sample,
                                                list(range(intval[0], intval[1], stepSize)))
        # remove duplicates
        positions_to_sample = np.unique(np.sort(positions_to_sample))
        if debug:
            print("sampling increased to {} from {}".format(
                len(positions_to_sample),
                orig_len))

    # skip regions that are filtered out
    if filter_out_tree:
        try:
            out_match = filter_out_tree.findOverlaps(chrom, start, end)
        except KeyError:
            out_match = []

        if len(out_match) > 0:
            for intval in out_match:
                positions_to_sample = \
                    positions_to_sample[
                        (positions_to_sample < intval[0]) |
                        (positions_to_sample >= intval[1])]
    return positions_to_sample


def countReadsPerGC_wrapper(args):
    return countReadsPerGC_worker(*args)


def countReadsPerGC_worker(chromNameBam,
                           start, end, stepSize, regionSize,
                           chrNameBamToBit, verbose=False):
    """given a genome region defined by
    (start, end), the GC content is quantified for
    regions of size regionSize that are contiguous
    """

    chromNameBit = chrNameBamToBit[chromNameBam]
    tbit = twobit.TwoBitFile(global_vars['2bit'])
    bam = bamHandler.openBam(global_vars['bam'])
    c = 1
    sub_reads_per_gc = []
    positions_to_sample = getPositionsToSample(chromNameBit,
                                               start, end, stepSize)

    for index in range(len(positions_to_sample)):
        i = positions_to_sample[index]
        # stop if region extends over the chromosome end
        if tbit.sequence_sizes()[chromNameBit] < i + regionSize:
            break

        try:
            gc = getGC_content(tbit[chromNameBit][i:i + regionSize])
        except Exception as detail:
            if verbose:
                print("{}:{}-{}".format(chromNameBit, i, i + regionSize))
                print(detail)
            continue
        numberReads = bam.count(chromNameBam, i, i + regionSize)
        sub_reads_per_gc.append((numberReads, gc))
        c += 1

    return sub_reads_per_gc


def tabulateGCcontent_wrapper(args):
    return tabulateGCcontent_worker(*args)


def tabulateGCcontent_worker(chromNameBam, start, end, stepSize,
                             fragmentLength,
                             chrNameBamToBit, verbose=False):
    r""" given genome regions, the GC content of the genome is tabulated for
    fragments of length 'fragmentLength' each 'stepSize' positions.

    >>> test = Tester()
    >>> args = test.testTabulateGCcontentWorker()
    >>> N_gc, F_gc = tabulateGCcontent_worker(*args)

    The forward read positions are:
    [1,  4,  10, 10, 16, 18]
    which correspond to a GC of
    [1,  1,  1,  1,  2,  1]

    The evaluated position are
    [0,  2,  4,  6,  8, 10, 12, 14, 16, 18]
    the corresponding GC is
    [2,  1,  1,  2,  2,  1,  2,  3,  2,  1]

    >>> print(N_gc)
    [0 4 5 1]
    >>> print(F_gc)
    [0 4 1 0]
    >>> test.set_filter_out_file()
    >>> chrNameBam2bit =  {'2L': 'chr2L'}

    Test for the filter out option
    >>> N_gc, F_gc = tabulateGCcontent_worker('2L', 0, 20, 2,
    ... {'median': 3}, chrNameBam2bit)
    >>> test.unset_filter_out_file()

    The evaluated positions are
    [ 0  2  8 10 12 14 16 18]
    >>> print(N_gc)
    [0 3 4 1]
    >>> print(F_gc)
    [0 3 1 0]

    Test for extra_sampling option
    >>> test.set_extra_sampling_file()
    >>> chrNameBam2bit =  {'2L': 'chr2L'}
    >>> res = tabulateGCcontent_worker('2L', 0, 20, 2,
    ... {'median': 3}, chrNameBam2bit)

    The new positions evaluated are
    [0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18]
    and the GC is
    [2, 1, 1, 0, 1, 2, 2, 1,  2,  3,  2,  1]
    >>> print(res[0])
    [1 5 5 1]
    >>> print(res[1])
    [0 5 1 0]

    """
    if start > end:
        raise NameError("start %d bigger that end %d" % (start, end))

    chromNameBit = chrNameBamToBit[chromNameBam]

    # array to keep track of the GC from regions of length 'fragmentLength'
    # from the genome. The index of the array is used to
    # indicate the gc content. The values inside the
    # array are counts. Thus, if N_gc[10] = 3, that means
    # that 3 regions have a gc_content of 10.
    subN_gc = np.zeros(fragmentLength['median'] + 1, dtype='int')
    subF_gc = np.zeros(fragmentLength['median'] + 1, dtype='int')

    tbit = twobit.TwoBitFile(global_vars['2bit'])
    bam = bamHandler.openBam(global_vars['bam'])
    peak = 0
    startTime = time.time()

    if verbose:
        print("[{:.3f}] computing positions to "
              "sample".format(time.time() - startTime))

    positions_to_sample = getPositionsToSample(chromNameBit,
                                               start, end, stepSize)

    read_counts = []
    # Optimize IO.
    # if the sample regions are far apart from each
    # other is faster to go to each location and fetch
    # the reads found there.
    # Otherwise, if the regions to sample are close to
    # each other, is faster to load all the reads in
    # a large region into memory and consider only
    # those falling into the positions to sample.
    # The following code gets the reads
    # that are at sampling positions that lie close together
    if np.mean(np.diff(positions_to_sample)) < 1000:
        start_pos = min(positions_to_sample)
        end_pos = max(positions_to_sample)
        if verbose:
            print("[{:.3f}] caching reads".format(time.time() - startTime))

        counts = np.bincount([r.pos - start_pos
                              for r in bam.fetch(chromNameBam, start_pos,
                                                 end_pos + 1)
                              if not r.is_reverse and r.pos >= start_pos],
                             minlength=end_pos - start_pos + 2)

        read_counts = counts[positions_to_sample - min(positions_to_sample)]
        if verbose:
            print("[{:.3f}] finish caching reads.".format(
                time.time() - startTime))

    countTime = time.time()

    c = 1
    for index in range(len(positions_to_sample)):
        i = positions_to_sample[index]
        # stop if the end of the chromosome is reached
        if i + fragmentLength['median'] > tbit.sequence_sizes()[chromNameBit]:
            break

        try:
            gc = getGC_content(
                tbit[chromNameBit][i:i + fragmentLength['median']],
                as_fraction=False)
        except Exception as detail:
            if verbose:
                print(detail)
            continue

        subN_gc[gc] += 1

        # count all reads at position 'i'
        if len(read_counts) == 0:  # case when no cache was done
            num_reads = len([x.pos for x in bam.fetch(chromNameBam, i, i + 1)
                             if x.is_reverse is False and x.pos == i])
        else:
            num_reads = read_counts[index]

        if num_reads >= global_vars['max_reads']:
            peak += 1
            continue

        subF_gc[gc] += num_reads
        if verbose:
            if index % 50000 == 0:
                endTime = time.time()
                print("%s processing %d (%.1f per sec) @ %s:%s-%s %s" %
                      (multiprocessing.current_process().name,
                       index, index / (endTime - countTime),
                       chromNameBit, start, end, stepSize))
        c += 1

    if verbose:
        endTime = time.time()
        print("%s processing %d (%.1f per sec) @ %s:%s-%s %s" %
              (multiprocessing.current_process().name,
               index, index / (endTime - countTime),
               chromNameBit, start, end, stepSize))
        print("%s total time %.1f @ %s:%s-%s %s" % (multiprocessing.current_process().name,
                                                    (endTime - startTime), chromNameBit, start, end, stepSize))

    return(subN_gc, subF_gc)


def tabulateGCcontent(fragmentLength, chrNameBitToBam, stepSize,
                      chromSizes, numberOfProcessors=None, verbose=False,
                      region=None):
    r"""
    Subdivides the genome or the reads into chunks to be analyzed in parallel
    using several processors. This codes handles the creation of
    workers that tabulate the GC content for small regions and then
    collects and integrates the results
    >>> test = Tester()
    >>> arg = test.testTabulateGCcontent()
    >>> res = tabulateGCcontent(*arg)
    >>> res
    array([[   0.        ,   18.        ,    1.        ],
           [   3.        ,   63.        ,    0.42857143],
           [   7.        ,  159.        ,    0.39622642],
           [  25.        ,  192.        ,    1.171875  ],
           [  28.        ,  215.        ,    1.17209302],
           [  16.        ,  214.        ,    0.6728972 ],
           [  12.        ,   95.        ,    1.13684211],
           [   9.        ,   24.        ,    3.375     ],
           [   3.        ,   11.        ,    2.45454545],
           [   0.        ,    0.        ,    1.        ],
           [   0.        ,    0.        ,    1.        ]])
    """
    global global_vars

    chrNameBamToBit = dict([(v, k) for k, v in chrNameBitToBam.items()])
    chunkSize = int(min(2e6, 4e5 / global_vars['reads_per_bp']))
    chromSizes = [(k, v) for k, v in chromSizes if k in list(chrNameBamToBit.keys())]

    imap_res = mapReduce.mapReduce((stepSize,
                                    fragmentLength, chrNameBamToBit,
                                    verbose),
                                   tabulateGCcontent_wrapper,
                                   chromSizes,
                                   genomeChunkLength=chunkSize,
                                   numberOfProcessors=numberOfProcessors,
                                   region=region)

    for subN_gc, subF_gc in imap_res:
        try:
            F_gc += subF_gc
            N_gc += subN_gc
        except NameError:
            F_gc = subF_gc
            N_gc = subN_gc

    scaling = sum(N_gc) // sum(F_gc)

    R_gc = np.array([float(F_gc[x]) / N_gc[x] * scaling
                     if N_gc[x] and F_gc[x] > 0 else 1
                     for x in range(len(F_gc))])

    data = np.transpose(np.vstack((F_gc, N_gc, R_gc)))
    return data


def countReadsPerGC(regionSize, chrNameBitToBam, stepSize,
                    chromSizes, numberOfProcessors=None, verbose=False,
                    region=None):
    r"""
    Computes for a region of size regionSize, the GC of the region
    and the number of reads that overlap it.
    >>> test = Tester()
    >>> arg = test.testCountReadsPerGC()
    >>> reads_per_gc = countReadsPerGC(*arg)
    >>> reads_per_gc[0:5,:]
    array([[ 132.        ,    0.44      ],
           [ 132.        ,    0.44      ],
           [ 133.        ,    0.44      ],
           [ 134.        ,    0.43666667],
           [ 134.        ,    0.44      ]])
    """
    global global_vars

    chrNameBamToBit = dict([(v, k) for k, v in chrNameBitToBam.items()])
    chunkSize = int(min(2e6, 4e5 / global_vars['reads_per_bp']))

    imap_res = mapReduce.mapReduce((stepSize,
                                    regionSize, chrNameBamToBit,
                                    verbose),
                                   countReadsPerGC_wrapper,
                                   chromSizes,
                                   genomeChunkLength=chunkSize,
                                   numberOfProcessors=numberOfProcessors,
                                   region=region)

    reads_per_gc = []
    for sub_reads_per_gc in imap_res:
        reads_per_gc += sub_reads_per_gc

    reads_per_gc = np.asarray(reads_per_gc)
    return reads_per_gc


def smooth(x, window_len=3):
    """
    *CURRENTLY* not being used
    smooths the values from the frequencies by taking the average
    of 'window_len' values.  window_len has to be an odd number
    """
    # do not smooth small arrays
    if len(x) < window_len * 2:
        return x
    i = 0
    y = x[:]
    half_width = (window_len - 1) / 2
    for i in range(0, len(x)):
        if i < half_width or i + half_width + 1 > len(x):
            continue
        else:
            y[i] = np.mean(x[i - half_width:i + half_width + 1])
    # clip low values, this avoid problems with zeros
    return y


def bin_by(x, y, nbins=10):
    """
    Bin x by y.
    Returns the binned "x" values and the left edges of the bins
    """
    bins = np.linspace(0, 1, nbins + 1)
    # To avoid extra bin for the max value
    bins[-1] += 1

    indices = np.digitize(y, bins)

    output = []
    for i in range(1, len(bins)):
        output.append(x[indices == i])

    # Just return the left edges of the bins
    bins = bins[:-1]

    return output, bins


def plotGCbias(file_name, frequencies, reads_per_gc, region_size, image_format=None):
    from matplotlib import use
    use('Agg')
    import matplotlib.pyplot as plt

    # prepare data for boxplot
    reads, GC = reads_per_gc.T
    reads_per_gc, bin_labels = bin_by(reads, GC, nbins=100)
    to_keep = [idx for idx, x in enumerate(bin_labels) if 0.2 <= x <= 0.7]
    reads_per_gc = [reads_per_gc[x] for x in to_keep]
    bin_labels = [bin_labels[x] for x in to_keep]

    title = "reads per regions of {} bp".format(region_size)
    fig = plt.figure(figsize=(6, 8))
    ax1 = fig.add_subplot(211, title=title)
    ax2 = fig.add_subplot(212,
                          title='normalized observed/expected read counts')

    # make boxplot

    bp = ax1.boxplot(reads_per_gc, notch=0, patch_artist=True)
    plt.setp(bp['boxes'], color='black', facecolor='LightGreen')
    plt.setp(bp['medians'], color='black')
    plt.setp(bp['whiskers'], color='black', linestyle='dashed')
    plt.setp(bp['fliers'], marker='None')
    # get the whisker that spands the most
    y_max = max([x.get_data()[1][1] for x in bp['whiskers']])
    ax1.set_ylim(0 - (y_max * 0.05), y_max * 1.05)
    ax1.set_ylabel('Number of reads')
    ax1.set_xlabel('GC fraction')

    xticks = [idx for idx, x in enumerate(bin_labels) if int(x * 100) % 10 == 0]

    ax1.set_xticks(xticks)
    ax1.set_xticklabels(["{:.1f}".format(bin_labels[x]) for x in xticks])

    x = np.linspace(0, 1, frequencies.shape[0])
    ax2.plot(x, np.log2(frequencies[:, 2]), color='#8c96f0')
    ax2.set_xlabel('GC fraction')
    ax2.set_ylabel('log2ratio observed/expected')
    ax2.set_xlim(0.2, 0.7)
    plt.tight_layout()
    plt.savefig(file_name, bbox_inches='tight', dpi=100, format=image_format)
    plt.close()


def main(args=None):
    args = parse_arguments().parse_args(args)

    if args.extraSampling:
        extra_sampling_file = args.extraSampling.name
        args.extraSampling.close()
    else:
        extra_sampling_file = None

    global global_vars
    global_vars = {}
    global_vars['2bit'] = args.genome
    global_vars['bam'] = args.bamfile
    global_vars['filter_out'] = args.blackListFileName
    global_vars['extra_sampling_file'] = extra_sampling_file

    tbit = twobit.TwoBitFile(global_vars['2bit'])
    bam = bamHandler.openBam(global_vars['bam'])

    if args.fragmentLength:
        fragment_len_dict = \
            {'median': args.fragmentLength}

    else:
        fragment_len_dict, __ = \
            get_read_and_fragment_length(args.bamfile, None,
                                         numberOfProcessors=args.numberOfProcessors,
                                         verbose=args.verbose)
        if not fragment_len_dict:
            print("\nPlease provide the fragment length used for the "
                  "sample preparation.\n")
            exit(1)

        fragment_len_dict = {'median': int(fragment_len_dict['median'])}

    chrNameBitToBam = tbitToBamChrName(list(tbit.sequence_sizes().keys()), bam.references)

    global_vars['genome_size'] = sum(tbit.sequence_sizes().values())
    global_vars['total_reads'] = bam.mapped
    global_vars['reads_per_bp'] = \
        float(global_vars['total_reads']) / args.effectiveGenomeSize

    confidence_p_value = float(1) / args.sampleSize

    # chromSizes: list of tuples
    chromSizes = [(bam.references[i], bam.lengths[i])
                  for i in range(len(bam.references))]

    # use poisson distribution to identify peaks that should be discarted.
    # I multiply by 4, because the real distribution of reads
    # vary depending on the gc content
    # and the global number of reads per bp may a be too low.
    # empirically, a value of at least 4 times as big as the
    # reads_per_bp was found.
    # Similarly for the min value, I divide by 4.
    global_vars['max_reads'] = \
        poisson(4 * global_vars['reads_per_bp'] *
                fragment_len_dict['median']).isf(confidence_p_value)
    # this may be of not use, unless the depth of sequencing is really high
    # as this value is close to 0
    global_vars['min_reads'] = \
        poisson(0.25 * global_vars['reads_per_bp'] *
                fragment_len_dict['median']).ppf(confidence_p_value)

    for key in global_vars:
        print("{}: {}".format(key, global_vars[key]))

    print("computing frequencies")
    # the GC of the genome is sampled each stepSize bp.
    stepSize = max(int(global_vars['genome_size'] / args.sampleSize), 1)
    print("stepSize: {}".format(stepSize))
    data = tabulateGCcontent(fragment_len_dict,
                             chrNameBitToBam, stepSize,
                             chromSizes,
                             numberOfProcessors=args.numberOfProcessors,
                             verbose=args.verbose,
                             region=args.region)

    np.savetxt(args.GCbiasFrequenciesFile.name, data)

    if args.biasPlot:
        reads_per_gc = countReadsPerGC(args.regionSize,
                                       chrNameBitToBam, stepSize * 10,
                                       chromSizes,
                                       numberOfProcessors=args.numberOfProcessors,
                                       verbose=args.verbose,
                                       region=args.region)
        plotGCbias(args.biasPlot, data, reads_per_gc, args.regionSize, image_format=args.plotFileFormat)


class Tester():
    def __init__(self):
        import os
        self.root = os.path.dirname(os.path.abspath(__file__)) + "/test/test_corrGC/"
        self.tbitFile = self.root + "sequence.2bit"
        self.bamFile = self.root + "test.bam"
        self.mappability = self.root + "mappability.bw"
        self.chrNameBam = '2L'
        self.chrNameBit = 'chr2L'
        bam = bamHandler.openBam(self.bamFile)
        tbit = twobit.TwoBitFile(self.tbitFile)
        global debug
        debug = 0
        global global_vars
        global_vars = {'2bit': self.tbitFile,
                       'bam': self.bamFile,
                       'filter_out': None,
                       'mappability': self.mappability,
                       'extra_sampling_file': None,
                       'max_reads': 5,
                       'min_reads': 0,
                       'min_reads': 0,
                       'reads_per_bp': 0.3,
                       'total_reads': bam.mapped,
                       'genome_size': sum(tbit.sequence_sizes().values())
                       }

    def testTabulateGCcontentWorker(self):
        stepSize = 2
        fragmentLength = {'min': 1, 'median': 3, 'max': 5}
        start = 0
        end = 20
        chrNameBam2bit = {'2L': 'chr2L'}
        return (self.chrNameBam,
                start, end, stepSize, fragmentLength, chrNameBam2bit)

    def set_filter_out_file(self):
        global global_vars
        global_vars['filter_out'] = self.root + "filter_out.bed"

    def unset_filter_out_file(self):
        global global_vars
        global_vars['filter_out'] = None

    def set_extra_sampling_file(self):
        global global_vars
        global_vars['extra_sampling_file'] = self.root + "extra_sampling.bed"

    def testTabulateGCcontent(self):
        fragmentLength = {'median': 10}
        chrNameBitToBam = {'chr2L': '2L'}
        stepSize = 1
        bam = bamHandler.openBam(global_vars['bam'])
        chromSizes = [(bam.references[i], bam.lengths[i])
                      for i in range(len(bam.references))]
        return (fragmentLength,
                chrNameBitToBam, stepSize, chromSizes, 1)

    def testCountReadsPerGC(self):
        regionSize = 300
        chrNameBitToBam = {'chr2L': '2L'}
        stepSize = 1
        bam = bamHandler.openBam(global_vars['bam'])
        chromSizes = [(bam.references[i], bam.lengths[i])
                      for i in range(len(bam.references))]
        return (regionSize,
                chrNameBitToBam, stepSize, chromSizes, 1)

if __name__ == "__main__":
    main()
