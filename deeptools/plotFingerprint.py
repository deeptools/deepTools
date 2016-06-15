#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
from matplotlib import use as mplt_use
mplt_use('Agg')
import matplotlib.pyplot as plt

import deeptools.countReadsPerBin as countR
from deeptools import parserCommon

old_settings = np.seterr(all='ignore')


def parse_arguments(args=None):
    parent_parser = parserCommon.getParentArgParse(binSize=False)
    required_args = get_required_args()
    output_args = get_output_args()
    optional_args = get_optional_args()
    read_options_parser = parserCommon.read_options()
    parser = argparse.ArgumentParser(
        parents=[required_args, output_args, read_options_parser,
                 optional_args, parent_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='This tool samples indexed BAM files '
        'and plots a profile of cumulative read coverages for each. '
        'All reads overlapping a window (bin) of the '
        'specified length are counted; '
        'these counts are sorted '
        'and the cumulative sum is finally plotted. ',
        conflict_handler='resolve',
        usage='An example usage is: plotFingerprint -b treatment.bam control.bam '
        '-plot fingerprint.png',
        add_help=False)

    return parser


def process_args(args=None):

    args = parse_arguments().parse_args(args)

    if args.labels and len(args.bamfiles) != len(args.labels):
        print("The number of labels does not match the number of BAM files.")
        exit(0)

    if not args.labels:
        args.labels = args.bamfiles

    return args


def get_required_args():
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    # define the arguments
    required.add_argument('--bamfiles', '-b',
                          metavar='bam files',
                          nargs='+',
                          help='List of indexed BAM files',
                          required=True)
    return parser


def get_optional_args():
    parser = argparse.ArgumentParser(add_help=False,
                                     conflict_handler='resolve')
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")

    optional.add_argument('--labels', '-l',
                          metavar='',
                          help='List of labels to use in the output. '
                          'If not given, the file names will be used instead. '
                          'Separate the labels by spaces.',
                          nargs='+')

    optional.add_argument('--binSize', '-bs',
                          help='Window size in base pairs to '
                          'sample the genome.',
                          default=500,
                          type=int)

    optional.add_argument('--numberOfSamples', '-n',
                          help='Number of bins that sampled from the genome, '
                          'for which the overlapping number of reads is computed.',
                          default=5e5,
                          type=int)

    optional.add_argument('--plotFileFormat',
                          metavar='',
                          help='image format type. If given, this option '
                          'overrides the image format based on the ending '
                          'given via --plotFile '
                          'ending. The available options are: "png", '
                          '"eps", "pdf" and "svg"',
                          choices=['png', 'pdf', 'svg', 'eps'])

    optional.add_argument('--plotTitle', '-T',
                          help='Title of the plot, to be printed on top of '
                          'the generated image. Leave blank for no title.',
                          default='')

    optional.add_argument('--skipZeros',
                          help='If set, then regions with zero overlapping reads'
                          'for *all* given BAM files are ignored. This '
                          'will result in a reduced number of read '
                          'counts than that specified in --numberOfSamples',
                          action='store_true')

    optional.add_argument('--outQualityMetrics',
                          help='Quality metrics can optionally be output to '
                          'this file. The file will have one row per input BAM '
                          'file and columns containing the following (in '
                          'order): area under the curve, X-intersect, and elbow '
                          'position (maximum of the second derivative)',
                          metavar='FILE.txt',
                          type=argparse.FileType('w'))

    optional.add_argument('--KLDsample',
                          help='Reference sample against which to compute the '
                          'Kullback-Leibler divergence. If this is not specified, '
                          'the divergence will not be calculated. If '
                          '--outQualityMetrics is not specified then this will '
                          'be ignored. Note that the KL divergence is calculated '
                          'on the log10 distribution of coverage counts, grouped '
                          'into 100 equally spaced bins. The input distribution '
                          'is modified slightly such that the resulting PMF '
                          'never goes to zero.',
                          metavar='sample.bam')

    return parser


def get_output_args():
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Output')
    group.add_argument('--plotFile', '-plot',
                       help='File name of the output figure. The file '
                       'ending will be used to determine the image '
                       'format. The available options are typically: "png", '
                       '"eps", "pdf" and "svg", e.g. : fingerprint.png.',
                       metavar='',
                       type=argparse.FileType('w'),
                       required=True)

    group.add_argument('--outRawCounts',
                       help='Output file name to save the read counts per bin.',
                       metavar='',
                       type=argparse.FileType('w'))

    return parser


def getKLD(args, idx, mat):
    """
    Computes the Kullback-Leibler divergence between two samples. The divergence
    between a sample and itself is "NA". If the reference sample doesn't exist
    then "NA" is returned.

    args: The input arguments
    idx:  The column index of the current sample
    mat:  The matrix of counts
    """

    # Get the index of the reference sample
    if args.KLDsample not in args.bamfiles:
        return "NA"
    refIdx = args.bamfiles.index(args.KLDsample)
    if refIdx == idx:
        return "NA"

    # Generate PMFs
    refVals = np.log10(mat[:, refIdx] + 0.1)
    sampleVals = np.log10(mat[:, idx] + 0.1)
    samplePMF, _ = np.histogram(sampleVals, bins=100, density=True)
    refPMF, _ = np.histogram(refVals, bins=100, density=True)

    # The reference PMF needs to be offset so q is never 0
    # The integral should still equal 1
    binSize = 1.0 / np.sum(refPMF)
    refPMF += 0.001
    integral = np.sum(refPMF) * binSize
    refPMF /= integral

    # Compute the KL divergence
    divergence = 0.0
    for i, (q, p) in enumerate(zip(refPMF, samplePMF)):
        if p > 0.0 and q > 0.0:
            divergence += p * np.log(p / q)

    return divergence


def main(args=None):
    args = process_args(args)

    cr = countR.CountReadsPerBin(
        args.bamfiles,
        args.binSize,
        args.numberOfSamples,
        blackListFileName=args.blackListFileName,
        numberOfProcessors=args.numberOfProcessors,
        verbose=args.verbose,
        region=args.region,
        extendReads=args.extendReads,
        minMappingQuality=args.minMappingQuality,
        ignoreDuplicates=args.ignoreDuplicates,
        center_read=args.centerReads,
        samFlag_include=args.samFlagInclude,
        samFlag_exclude=args.samFlagExclude)

    num_reads_per_bin = cr.run()
    if num_reads_per_bin.sum() == 0:
        import sys
        sys.stderr.write(
            "\nNo reads were found in {} regions sampled. Check that the\n"
            "min mapping quality is not overly high and that the \n"
            "chromosome names between bam files are consistant.\n"
            "\n".format(num_reads_per_bin.shape[0]))
        exit(1)

    if args.skipZeros:
        num_reads_per_bin = countR.remove_row_of_zeros(num_reads_per_bin)

    total = len(num_reads_per_bin[:, 0])
    x = np.arange(total).astype('float') / total  # normalize from 0 to 1

    i = 0
    # matplotlib won't iterate through line styles by itself
    pyplot_line_styles = sum([7 * ["-"], 7 * ["--"], 7 * ["-."], 7 * [":"], 7 * ["."]], [])
    for i, reads in enumerate(num_reads_per_bin.T):
        count = np.cumsum(np.sort(reads))
        count = count / count[-1]  # to normalize y from 0 to 1
        j = i % 35
        plt.plot(x, count, label=args.labels[i], linestyle=pyplot_line_styles[j])
        plt.xlabel('rank')
        plt.ylabel('fraction w.r.t. bin with highest coverage')
    plt.legend(loc='upper left')
    plt.suptitle(args.plotTitle)
    # set the plotFileFormat explicitly to None to trigger the
    # format from the file-extension
    if not args.plotFileFormat:
        args.plotFileFormat = None

    plt.savefig(args.plotFile.name, bbox_inches=0, format=args.plotFileFormat)
    plt.close()

    if args.outRawCounts:
        args.outRawCounts.write("'" + "'\t'".join(args.labels) + "'\n")
        fmt = "\t".join(np.repeat('%d', num_reads_per_bin.shape[1])) + "\n"
        for row in num_reads_per_bin:
            args.outRawCounts.write(fmt % tuple(row))
        args.outRawCounts.close()

    if args.outQualityMetrics:
        args.outQualityMetrics.write("Sample\tAUC\tX-intercept\tElbow Point")
        if args.KLDsample:
            args.outQualityMetrics.write("\tKL Divergence")
        args.outQualityMetrics.write("\n")
        line = np.arange(num_reads_per_bin.shape[0]) / float(num_reads_per_bin.shape[0] - 1)
        for idx, reads in enumerate(num_reads_per_bin.T):
            counts = np.cumsum(np.sort(reads))
            counts = counts / float(counts[-1])
            AUC = np.sum(counts)
            XInt = (np.argmax(counts > 0) + 1) / float(counts.shape[0])
            elbow = (np.argmax(line - counts) + 1) / float(counts.shape[0])
            if args.KLDsample:
                KLD = getKLD(args, idx, num_reads_per_bin)
                args.outQualityMetrics.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(args.labels[idx], AUC, XInt, elbow, KLD))
            else:
                args.outQualityMetrics.write("{0}\t{1}\t{2}\t{3}\n".format(args.labels[idx], AUC, XInt, elbow))
        args.outQualityMetrics.close()

if __name__ == "__main__":
    main()
