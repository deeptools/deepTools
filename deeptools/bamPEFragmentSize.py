#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import numpy as np

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
from deeptools import cm  # noqa: F401
import matplotlib.pyplot as plt

import plotly.offline as py
import plotly.graph_objs as go

# own tools
from deeptools.parserCommon import writableFile
from deeptools.getFragmentAndReadSize import get_read_and_fragment_length
from deeptools._version import __version__


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='This tool calculates the fragment sizes for read pairs given a BAM file from paired-end sequencing.'
        'Several regions are sampled depending on the '
        'size of the genome and number of processors to estimate the'
        'summary statistics on the fragment lengths. '
        'Properly paired reads are preferred for computation, i.e., '
        'it will only use discordant pairs if no concordant alignments '
        'overlap with a given region. '
        'The default setting simply prints the summary statistics to the screen.')
    parser.add_argument('--bamfiles', '-b',
                        help='List of BAM files to process',
                        nargs='+',
                        metavar='bam files')

    parser.add_argument('--histogram', '-hist', '-o',
                        help='Save a .png file with a histogram '
                        'of the fragment length distribution.',
                        metavar='FILE')

    parser.add_argument('--plotFileFormat',
                        metavar='FILETYPE',
                        help='Image format type. If given, this option '
                        'overrides the image format based on the plotFile '
                        'ending. The available options are: png, '
                        'eps, pdf, svg and plotly.',
                        default=None,
                        choices=['png', 'pdf', 'svg', 'eps', 'plotly'])

    parser.add_argument('--numberOfProcessors', '-p',
                        help='Number of processors to use. The default is '
                        'to use 1. (Default: %(default)s)',
                        metavar="INT",
                        type=int,
                        default=1,
                        required=False)
    parser.add_argument('--samplesLabel',
                        help='Labels for the samples plotted. The '
                        'default is to use the file name of the '
                        'sample. The sample labels should be separated '
                        'by spaces and quoted if a label itself'
                        'contains a space E.g. --samplesLabel label-1 "label 2"  ',
                        nargs='+')
    parser.add_argument('--plotTitle', '-T',
                        help='Title of the plot, to be printed on top of '
                        'the generated image. Leave blank for no title. (Default: %(default)s)',
                        default='')
    parser.add_argument('--maxFragmentLength',
                        help='The maximum fragment length in the histogram. A value of 0 (the default) indicates to use twice the mean fragment length. (Default: %(default)s)',
                        default=0,
                        type=int)
    parser.add_argument('--logScale',
                        help='Plot on the log scale',
                        action='store_true')
    parser.add_argument('--binSize', '-bs',
                        metavar='INT',
                        help='Length in bases of the window used to sample the genome. (Default: %(default)s)',
                        default=1000,
                        type=int)
    parser.add_argument('--distanceBetweenBins', '-n',
                        metavar='INT',
                        help='To reduce the computation time, not every possible genomic '
                        'bin is sampled. This option allows you to set the distance '
                        'between bins actually sampled from. Larger numbers are sufficient '
                        'for high coverage samples, while smaller values are useful for '
                        'lower coverage samples. Note that if you specify a value that '
                        'results in too few (<1000) reads sampled, the value will be '
                        'decreased. (Default: %(default)s)',
                        default=1000000,
                        type=int)
    parser.add_argument('--blackListFileName', '-bl',
                        help="A BED file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered.",
                        metavar="BED file",
                        required=False)
    parser.add_argument('--table',
                        metavar='FILE',
                        help='In addition to printing read and fragment length metrics to the screen, write them to the given file in tabular format.',
                        required=False)
    parser.add_argument('--outRawFragmentLengths',
                        metavar='FILE',
                        required=False,
                        type=writableFile,
                        help='Save the fragment (or read if the input is single-end) length and their associated number of occurrences to a tab-separated file. Columns are length, number of occurrences, and the sample label.')
    parser.add_argument('--verbose',
                        help='Set if processing data messages are wanted.',
                        action='store_true',
                        required=False)
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def getDensity(lengths, minVal, maxVal):
    """
    This is essentially computing what hist() in matplotlib is doing and returning the results.
    This then allows us to free up the memory consumed by each sample rather than returning it all back to main() for plotting.
    """
    n, bins, patches = plt.hist(lengths, bins=100, range=(minVal, maxVal), density=True)
    plt.clf()
    return (n, bins)


def getFragSize(bam, args, idx, outRawFrags):
    fragment_len_dict, read_len_dict = get_read_and_fragment_length(bam, return_lengths=True,
                                                                    blackListFileName=args.blackListFileName,
                                                                    numberOfProcessors=args.numberOfProcessors,
                                                                    verbose=args.verbose,
                                                                    binSize=args.binSize,
                                                                    distanceBetweenBins=args.distanceBetweenBins)

    if outRawFrags:
        label = bam
        if args.samplesLabel and idx < len(args.samplesLabel):
            label = args.samplesLabel[idx]
        if fragment_len_dict:
            fragment_len_dict['lengths'] = [int(x) for x in fragment_len_dict['lengths']]
            cnts = np.bincount(fragment_len_dict['lengths'], minlength=int(fragment_len_dict['max']) + 1)
        else:
            read_len_dict['lengths'] = [int(x) for x in read_len_dict['lengths']]
            cnts = np.bincount(read_len_dict['lengths'], minlength=int(read_len_dict['max']) + 1)
        for idx, v in enumerate(cnts):
            if v > 0:
                outRawFrags.write("{}\t{}\t{}\n".format(idx, v, label))

    if args.samplesLabel and idx < len(args.samplesLabel):
        print("\n\nSample label: {}".format(args.samplesLabel[idx]))
    else:
        print("\n\nBAM file : {}".format(bam))

    if fragment_len_dict:
        if fragment_len_dict['mean'] == 0:
            print("No pairs were found. Is the data from a paired-end sequencing experiment?")

        print("Sample size: {}\n".format(fragment_len_dict['sample_size']))

        print("Fragment lengths:")
        print("Min.: {}\n1st Qu.: {}\nMean: {}\nMedian: {}\n"
              "3rd Qu.: {}\nMax.: {}\nStd: {}".format(fragment_len_dict['min'],
                                                      fragment_len_dict['qtile25'],
                                                      fragment_len_dict['mean'],
                                                      fragment_len_dict['median'],
                                                      fragment_len_dict['qtile75'],
                                                      fragment_len_dict['max'],
                                                      fragment_len_dict['std']))
        print("MAD: {}\nLen. 10%: {}\nLen. 20%: {}\nLen. 30%: {}\nLen. 40%: {}\nLen. 60%: {}\nLen. 70%: {}\nLen. 80%: {}\nLen. 90%: {}\nLen. 99%: {}\n".format(fragment_len_dict['mad'],
                                                                                                                                                               fragment_len_dict['qtile10'],
                                                                                                                                                               fragment_len_dict['qtile20'],
                                                                                                                                                               fragment_len_dict['qtile30'],
                                                                                                                                                               fragment_len_dict['qtile40'],
                                                                                                                                                               fragment_len_dict['qtile60'],
                                                                                                                                                               fragment_len_dict['qtile70'],
                                                                                                                                                               fragment_len_dict['qtile80'],
                                                                                                                                                               fragment_len_dict['qtile90'],
                                                                                                                                                               fragment_len_dict['qtile99']))
    else:
        print("No pairs were found. Is the data from a paired-end sequencing experiment?")

    print("\nRead lengths:")
    print("Sample size: {}\n".format(read_len_dict['sample_size']))
    print("Min.: {}\n1st Qu.: {}\nMean: {}\nMedian: {}\n"
          "3rd Qu.: {}\nMax.: {}\nStd: {}".format(read_len_dict['min'],
                                                  read_len_dict['qtile25'],
                                                  read_len_dict['mean'],
                                                  read_len_dict['median'],
                                                  read_len_dict['qtile75'],
                                                  read_len_dict['max'],
                                                  read_len_dict['std']))
    print("MAD: {}\nLen. 10%: {}\nLen. 20%: {}\nLen. 30%: {}\nLen. 40%: {}\nLen. 60%: {}\nLen. 70%: {}\nLen. 80%: {}\nLen. 90%: {}\nLen. 99%: {}\n".format(read_len_dict['mad'],
                                                                                                                                                           read_len_dict['qtile10'],
                                                                                                                                                           read_len_dict['qtile20'],
                                                                                                                                                           read_len_dict['qtile30'],
                                                                                                                                                           read_len_dict['qtile40'],
                                                                                                                                                           read_len_dict['qtile60'],
                                                                                                                                                           read_len_dict['qtile70'],
                                                                                                                                                           read_len_dict['qtile80'],
                                                                                                                                                           read_len_dict['qtile90'],
                                                                                                                                                           read_len_dict['qtile99']))

    # The read and fragment lists will just eat up memory if not removed!
    if args.histogram:
        if fragment_len_dict:
            maxVal = fragment_len_dict['mean'] * 2
            minVal = fragment_len_dict['min']
        else:
            maxVal = read_len_dict['mean'] * 2
            minVal = read_len_dict['min']
        if args.maxFragmentLength > 0:
            maxVal = args.maxFragmentLength

        if fragment_len_dict:
            fragment_len_dict['lengths'] = getDensity(fragment_len_dict['lengths'], minVal, maxVal)
        if read_len_dict:
            read_len_dict['lengths'] = getDensity(read_len_dict['lengths'], minVal, maxVal)
    else:
        if fragment_len_dict:
            del fragment_len_dict['lengths']
        if read_len_dict:
            del read_len_dict['lengths']

    return (fragment_len_dict, read_len_dict)


def printTable(args, fragDict, readDict):
    """
    Print the read and fragment dictionary in more easily parsable tabular format to a file.
    """
    of = open(args.table, "w")
    of.write("\tFrag. Sampled")
    of.write("\tFrag. Len. Min.\tFrag. Len. 1st. Qu.\tFrag. Len. Mean\tFrag. Len. Median\tFrag. Len. 3rd Qu.\tFrag. Len. Max\tFrag. Len. Std.")
    of.write("\tFrag. Med. Abs. Dev.\tFrag. Len. 10%\tFrag. Len. 20%\tFrag. Len. 30%\tFrag. Len. 40%\tFrag. Len. 60%\tFrag. Len. 70%\tFrag. Len. 80%\tFrag. Len. 90%\tFrag. Len. 99%")
    of.write("\tReads Sampled")
    of.write("\tRead Len. Min.\tRead Len. 1st. Qu.\tRead Len. Mean\tRead Len. Median\tRead Len. 3rd Qu.\tRead Len. Max\tRead Len. Std.")
    of.write("\tRead Med. Abs. Dev.\tRead Len. 10%\tRead Len. 20%\tRead Len. 30%\tRead Len. 40%\tRead Len. 60%\tRead Len. 70%\tRead Len. 80%\tRead Len. 90%\tRead Len. 99%\n")

    for idx, bam in enumerate(args.bamfiles):
        if args.samplesLabel and idx < len(args.samplesLabel):
            of.write(args.samplesLabel[idx])
        else:
            of.write(bam)
        if fragDict is not None and fragDict[bam] is not None:
            d = fragDict[bam]
            of.write("\t{}".format(d['sample_size']))
            of.write("\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(d['min'],
                                                           d['qtile25'],
                                                           d['mean'],
                                                           d['median'],
                                                           d['qtile75'],
                                                           d['max'],
                                                           d['std']))
            of.write("\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(d['mad'],
                                                                       d['qtile10'],
                                                                       d['qtile20'],
                                                                       d['qtile30'],
                                                                       d['qtile40'],
                                                                       d['qtile60'],
                                                                       d['qtile70'],
                                                                       d['qtile80'],
                                                                       d['qtile90'],
                                                                       d['qtile99']))
        else:
            of.write("\t0")
            of.write("\t0\t0\t0\t0\t0\t0\t0")
            of.write("\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0")
        d = readDict[bam]
        of.write("\t{}".format(d['sample_size']))
        of.write("\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(d['min'],
                                                       d['qtile25'],
                                                       d['mean'],
                                                       d['median'],
                                                       d['qtile75'],
                                                       d['max'],
                                                       d['std']))
        of.write("\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(d['mad'],
                                                                     d['qtile10'],
                                                                     d['qtile20'],
                                                                     d['qtile30'],
                                                                     d['qtile40'],
                                                                     d['qtile60'],
                                                                     d['qtile70'],
                                                                     d['qtile80'],
                                                                     d['qtile90'],
                                                                     d['qtile99']))
    of.close()


def main(args=None):
    args = parse_arguments().parse_args(args)

    fraglengths = {}
    readlengths = {}
    of = None
    if args.outRawFragmentLengths is not None:
        of = open(args.outRawFragmentLengths, "w")
        of.write("#bamPEFragmentSize\nSize\tOccurrences\tSample\n")
    for idx, bam in enumerate(args.bamfiles):
        f, r = getFragSize(bam, args, idx, of)
        fraglengths[bam] = f
        readlengths[bam] = r

    if args.table is not None:
        printTable(args, fraglengths, readlengths)

    if args.histogram:
        if args.samplesLabel:
            if len(args.bamfiles) != len(args.samplesLabel):
                sys.exit("The number of labels does not match the number of BAM files.")
            else:
                labels = args.samplesLabel
        else:
            labels = list(fraglengths.keys())

        i = 0
        data = []
        for bam in fraglengths.keys():
            d = fraglengths[bam]
            if d is None:
                d = readlengths[bam]
            if args.maxFragmentLength > 0:
                maxVal = args.maxFragmentLength
            else:
                maxVal = d['mean'] * 2

            if args.plotFileFormat == 'plotly':
                trace = go.Histogram(x=d['lengths'],
                                     histnorm='probability',
                                     opacity=0.5,
                                     name=labels[i],
                                     nbinsx=100,
                                     xbins=dict(start=d['min'], end=maxVal))
                data.append(trace)
            else:
                plt.bar(d['lengths'][1][:-1], height=d['lengths'][0],
                        width=d['lengths'][1][1:] - d['lengths'][1][:-1],
                        align='edge', log=args.logScale,
                        alpha=0.5, label=labels[i])
            i += 1

        if args.plotFileFormat == 'plotly':
            fig = go.Figure()
            fig.add_traces(data)
            fig['layout']['yaxis1'].update(title='Frequency')
            fig['layout']['xaxis1'].update(title='Fragment Length')
            fig['layout'].update(title=args.plotTitle)
            fig['layout'].update(showlegend=True)
            if args.logScale:
                fig['layout']['yaxis1'].update(type='log')
            py.plot(fig, filename=args.histogram, auto_open=False)
        else:
            plt.xlabel('Fragment Length')
            plt.ylabel('Frequency')
            plt.legend(loc='upper right')
            plt.title(args.plotTitle)
            plt.savefig(args.histogram, bbox_inches=0, format=args.plotFileFormat)
            plt.close()


if __name__ == "__main__":
    main()
