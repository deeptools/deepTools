#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np

import deeptools.countReadsPerBin as countR
from deeptools import parserCommon
from deeptools._version import __version__

old_settings = np.seterr(all='ignore')


def parse_arguments(args=None):
    parser = \
        argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description="""

``multiBamSummary`` computes the read coverages for genomic regions for typically two or more BAM files.
The analysis can be performed for the entire genome by running the program in 'bins' mode.
If you want to count the read coverage for specific regions only, use the ``BED-file`` mode instead.
The standard output of ``multiBamSummary`` is a compressed numpy array (``.npz``).
It can be directly used to calculate and visualize pairwise correlation values between the read coverages using the tool 'plotCorrelation'.
Similarly, ``plotPCA`` can be used for principal component analysis of the read coverages using the .npz file.
Note that using a single bigWig file is only recommended if you want to produce a bedGraph file (i.e., with the ``--outRawCounts`` option; the default output file cannot be used by ANY deepTools program if only a single file was supplied!).

A detailed sub-commands help is available by typing:

  multiBamSummary bins -h

  multiBamSummary BED-file -h


""",
            epilog='example usages:\n'
                   'multiBamSummary bins --bamfiles file1.bam file2.bam -out results.npz \n\n'
                   'multiBamSummary BED-file --BED selection.bed --bamfiles file1.bam file2.bam \n'
                   '-out results.npz'
                   ' \n\n',
            conflict_handler='resolve')

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))
    subparsers = parser.add_subparsers(
        title="commands",
        dest='command',
        description='subcommands',
        help='subcommands',
        metavar='')

    parent_parser = parserCommon.getParentArgParse(binSize=False)
    read_options_parser = parserCommon.read_options()

    # bins mode options
    subparsers.add_parser(
        'bins',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[bamcorrelate_args(case='bins'),
                 parent_parser, read_options_parser,
                 parserCommon.gtf_options(suppress=True)
                 ],
        help="The coverage calculation is done for consecutive bins of equal "
             "size (10 kilobases by default). This mode is useful to assess the "
             "genome-wide similarity of BAM files. The bin size and "
             "distance between bins can be adjusted.",
        add_help=False,
        usage='%(prog)s '
              '--bamfiles file1.bam file2.bam '
              '-out results.npz \n')

    # BED file arguments
    subparsers.add_parser(
        'BED-file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[bamcorrelate_args(case='BED-file'),
                 parent_parser, read_options_parser,
                 parserCommon.gtf_options()
                 ],
        help="The user provides a BED file that contains all regions "
             "that should be considered for the coverage analysis. A "
             "common use is to compare ChIP-seq coverages between two "
             "different samples for a set of peak regions.",
        usage='%(prog)s --BED selection.bed --bamfiles file1.bam file2.bam -out results.npz\n',
        add_help=False)

    return parser


def bamcorrelate_args(case='bins'):
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    # define the arguments
    required.add_argument('--bamfiles', '-b',
                          metavar='FILE1 FILE2',
                          help='List of indexed bam files separated by spaces.',
                          nargs='+',
                          required=True)

    required.add_argument('--outFileName', '-out',
                          help='File name to save the coverage matrix. This matrix '
                               'can be subsequently plotted using plotCorrelation or '
                               'or plotPCA.',
                          required=True)

    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")
    optional.add_argument('--labels', '-l',
                          metavar='sample1 sample2',
                          help='User defined labels instead of default labels from '
                               'file names. '
                               'Multiple labels have to be separated by a space, e.g. '
                               '--labels sample1 sample2 sample3',
                          nargs='+')

    if case == 'bins':
        optional.add_argument('--binSize', '-bs',
                              metavar='INT',
                              help='Length in bases of the window used '
                                   'to sample the genome.',
                              default=10000,
                              type=int)

        optional.add_argument('--distanceBetweenBins', '-n',
                              metavar='INT',
                              help='By default, multiBamSummary considers consecutive '
                              'bins of the specified --binSize. However, to '
                              'reduce the computation time, a larger distance '
                              'between bins can by given. Larger distances '
                              'result in fewer bins considered.',
                              default=0,
                              type=int)

        required.add_argument('--BED',
                              help=argparse.SUPPRESS,
                              default=None)
    else:
        optional.add_argument('--binSize', '-bs',
                              help=argparse.SUPPRESS,
                              default=10000,
                              type=int)

        optional.add_argument('--distanceBetweenBins', '-n',
                              help=argparse.SUPPRESS,
                              metavar='INT',
                              default=0,
                              type=int)

        required.add_argument('--BED',
                              help='Limits the coverage analysis to '
                              'the regions specified in these files.',
                              metavar='FILE1.bed FILE2.bed',
                              nargs='+',
                              required=True)

    group = parser.add_argument_group('Output optional options')

    group.add_argument('--outRawCounts',
                       help='Save the counts per region to a tab-delimited file.',
                       metavar='FILE')

    return parser


def process_args(args=None):
    args = parse_arguments().parse_args(args)

    if args.labels and len(args.bamfiles) != len(args.labels):
        print("The number of does not match the number of bam files.")
        exit(0)
    if not args.labels:
        args.labels = [os.path.basename(x) for x in args.bamfiles]

    return args


def main(args=None):
    """
    1. get read counts at different positions either
    all of same length or from genomic regions from the BED file

    2. save data for further plotting

    """
    args = process_args(args)

    if 'BED' in args:
        bed_regions = args.BED
    else:
        bed_regions = None

    if len(args.bamfiles) == 1 and not args.outRawCounts:
        sys.stderr.write("You've input a single BAM file and not specified "
                         "--outRawCounts. The resulting output will NOT be "
                         "useful with any deepTools program!\n")

    stepsize = args.binSize + args.distanceBetweenBins
    c = countR.CountReadsPerBin(
        args.bamfiles,
        args.binSize,
        numberOfSamples=None,
        numberOfProcessors=args.numberOfProcessors,
        verbose=args.verbose,
        region=args.region,
        bedFile=bed_regions,
        blackListFileName=args.blackListFileName,
        extendReads=args.extendReads,
        minMappingQuality=args.minMappingQuality,
        ignoreDuplicates=args.ignoreDuplicates,
        center_read=args.centerReads,
        samFlag_include=args.samFlagInclude,
        samFlag_exclude=args.samFlagExclude,
        minFragmentLength=args.minFragmentLength,
        maxFragmentLength=args.maxFragmentLength,
        stepSize=stepsize,
        zerosToNans=False,
        out_file_for_raw_data=args.outRawCounts)

    num_reads_per_bin = c.run(allArgs=args)

    sys.stderr.write("Number of bins "
                     "found: {}\n".format(num_reads_per_bin.shape[0]))

    if num_reads_per_bin.shape[0] < 2:
        exit("ERROR: too few non zero bins found.\n"
             "If using --region please check that this "
             "region is covered by reads.\n")

    # numpy will append .npz to the file name if we don't do this...
    f = open(args.outFileName, "wb")
    np.savez_compressed(f,
                        matrix=num_reads_per_bin,
                        labels=args.labels)
    f.close()

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


if __name__ == "__main__":
    main()
