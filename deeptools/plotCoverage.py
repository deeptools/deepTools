#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np

from matplotlib import use as mplt_use
mplt_use('Agg')
import matplotlib.pyplot as plt

import deeptools.countReadsPerBin as countR
from deeptools import parserCommon
from deeptools._version import __version__

old_settings = np.seterr(all='ignore')


def parse_arguments(args=None):
    parent_parser = parserCommon.getParentArgParse(binSize=False)
    read_options_parser = parserCommon.read_options()

    parser = \
        argparse.ArgumentParser(
            parents=[required_args(), parent_parser, read_options_parser],
            formatter_class=argparse.RawDescriptionHelpFormatter,
            add_help=False,
            description="""

This tool is useful to assess the sequencing depth of a given sample.
It samples 1 million bp, counts the number of overlapping reads and can report
a histogram that tells you how many bases are covered how many times.
Multiple BAM files are accepted, but they all should correspond to the same genome assembly.

detailed usage help:
 $ plotCoverage  -h

""",
            epilog='example usages:\nplotCoverage '
                   '--bamfiles file1.bam file2.bam -out results.png\n\n'
                   ' \n\n',
            conflict_handler='resolve')

    parser.add_argument('--version', action='version',
                        version='plotCoverage {}'.format(__version__))

    return parser


def process_args(args=None):
    args = parse_arguments().parse_args(args)

    if args.labels and len(args.bamfiles) != len(args.labels):
        print("The number of labels does not match the number of BAM files.")
        exit(0)
    if not args.labels:
        args.labels = [os.path.basename(x) for x in args.bamfiles]

    return args


def required_args():
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    required.add_argument('--bamfiles', '-b',
                          metavar='FILE1 FILE2',
                          help='List of indexed BAM files separated by spaces.',
                          nargs='+',
                          required=True)

    required.add_argument('--plotFile', '-o',
                          help='File name to save the plot to.',
                          type=argparse.FileType('w'),
                          required=True)

    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")
    optional.add_argument('--labels', '-l',
                          metavar='sample1 sample2',
                          help='User defined labels instead of default labels from '
                               'file names. '
                               'Multiple labels have to be separated by spaces, e.g. '
                               '--labels sample1 sample2 sample3',
                          nargs='+')

    optional.add_argument('--plotTitle', '-T',
                          help='Title of the plot, to be printed on top of '
                          'the generated image. Leave blank for no title.',
                          default='')

    optional.add_argument('--skipZeros',
                          help='By setting this option, genomic regions '
                          'that have zero or nan values in _all_ samples '
                          'are excluded.',
                          action='store_true',
                          required=False)

    optional.add_argument('--numberOfSamples', '-n',
                          help='Number of 1 bp regions to sample. Default 1 million.',
                          required=False,
                          type=int,
                          default=1000000)

    optional.add_argument('--outRawCounts',
                          help='Save raw counts (coverages) to file.',
                          metavar='FILE')

    optional.add_argument('--plotFileFormat',
                          metavar='FILETYPE',
                          help='Image format type. If given, this option '
                          'overrides the image format based on the plotFile '
                          'ending. The available options are: png, '
                          'eps, pdf and svg.',
                          default=None,
                          choices=['png', 'pdf', 'svg', 'eps'])

    return parser


def main(args=None):
    args = process_args(args)
    cr = countR.CountReadsPerBin(args.bamfiles,
                                 binLength=1,
                                 numberOfSamples=args.numberOfSamples,
                                 numberOfProcessors=args.numberOfProcessors,
                                 verbose=args.verbose,
                                 region=args.region,
                                 blackListFileName=args.blackListFileName,
                                 extendReads=args.extendReads,
                                 minMappingQuality=args.minMappingQuality,
                                 ignoreDuplicates=args.ignoreDuplicates,
                                 center_read=args.centerReads,
                                 samFlag_include=args.samFlagInclude,
                                 samFlag_exclude=args.samFlagExclude,
                                 minFragmentLength=args.minFragmentLength,
                                 maxFragmentLength=args.maxFragmentLength,
                                 out_file_for_raw_data=args.outRawCounts)

    num_reads_per_bin = cr.run()

    sys.stderr.write("Number of non zero bins "
                     "used: {}\n".format(num_reads_per_bin.shape[0]))

    if args.outRawCounts:
        # append to the generated file the
        # labels
        header = "#'chr'\t'start'\t'end'\t"
        header += "'" + "'\t'".join(args.labels) + "'\n"
        f = open(args.outRawCounts, 'r+')
        content = f.read()
        f.seek(0, 0)
        f.write(header + content)
        f.close()

    if num_reads_per_bin.shape[0] < 2:
        exit("ERROR: too few non-zero bins found.\n"
             "If using --region please check that this "
             "region is covered by reads.\n")

    if args.skipZeros:
        num_reads_per_bin = countR.remove_row_of_zeros(num_reads_per_bin)

    fig, axs = plt.subplots(1, 2, figsize=(15, 5))
    plt.suptitle(args.plotTitle)
    # plot up to two std from mean
    num_reads_per_bin = num_reads_per_bin.astype(int)
    sample_mean = num_reads_per_bin.mean(axis=0)
    sample_std = num_reads_per_bin.std(axis=0)
    sample_max = num_reads_per_bin.max(axis=0)
    sample_min = num_reads_per_bin.min(axis=0)
    sample_25 = np.percentile(num_reads_per_bin, 25, axis=0)
    sample_50 = np.percentile(num_reads_per_bin, 50, axis=0)
    sample_75 = np.percentile(num_reads_per_bin, 75, axis=0)

    # use the largest 99th percentile from all samples to set the x_max value
    x_max = np.max(np.percentile(num_reads_per_bin, 99, axis=0))
    # plot coverage
    # print headers for text output
    print("sample\tmean\tstd\tmin\t25%\t50%\t75%\tmax")
    # the determination of a sensible value for y_max of the first plot (fraction of bases sampled vs.
    # coverage) is important because, depending on the data,
    # it becomes very difficult to see the lines in the plot. For example, if the coverage of a sample
    # is a nice gaussian curve with a large mean of 50. Then a sensible range for the y axis (fraction of
    # reads having coverage=x) is (0, 0.02) which nicely shows the coverage curve. If instead the coverage is
    # very por and centers close to 1 then a good y axis range is (0,1).

    # the current implementation aims to find the y_value for which 50% of the reads >= x (coverage) and
    # sets that as the x_axis range.
    y_max = []
    for idx, col in enumerate(num_reads_per_bin.T):
        frac_reads_per_coverage = np.bincount(col.astype(int)).astype(float) / num_reads_per_bin.shape[0]
        axs[0].plot(frac_reads_per_coverage, label="{}, mean={:.1f}".format(args.labels[idx], sample_mean[idx]))
        csum = np.bincount(col.astype(int))[::-1].cumsum()
        csum_frac = csum.astype(float)[::-1] / csum.max()
        axs[1].plot(csum_frac, label=args.labels[idx])
        # find the indexes (i.e. the x values) for which the cumulative distribution 'fraction of bases
        # sampled >= coverage' where fraction of bases sampled = 50%: `np.flatnonzero(csum_frac>0.5)`
        # then find the fraction of bases sampled that that have the largest x
        y_max.append(frac_reads_per_coverage[max(np.flatnonzero(csum_frac > 0.5))])
        print("{}\t{:0.2f}\t{:0.2f}\t{}\t{}\t{}\t{}\t{}\t".format(args.labels[idx],
                                                                  sample_mean[idx],
                                                                  sample_std[idx],
                                                                  sample_min[idx],
                                                                  sample_25[idx],
                                                                  sample_50[idx],
                                                                  sample_75[idx],
                                                                  sample_max[idx],
                                                                  ))

    # The 'good' x-axis is computed for each sample. The lower value is favored in which
    # distributions with a wider x-range can better be seen.
    y_max = min(y_max)
    axs[0].set_ylim(0, min(1, y_max + (y_max * 0.10)))
    axs[0].set_xlim(0, x_max)
    axs[0].set_xlabel('coverage (#reads per bp)')
    axs[0].legend(fancybox=True, framealpha=0.5)
    axs[0].set_ylabel('fraction of bases sampled')
    # plot cumulative coverage
    axs[1].set_xlim(0, x_max)
    axs[1].set_xlabel('coverage (#reads per bp)')
    axs[1].set_ylabel('fraction of bases sampled >= coverage')
    axs[1].legend(fancybox=True, framealpha=0.5)
    plt.savefig(args.plotFile.name, format=args.plotFileFormat)
    plt.close()

if __name__ == "__main__":
    main()
