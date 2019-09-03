#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import os.path
import numpy as np
import multiprocessing
from deeptools import parserCommon
from deeptools._version import __version__
from deeptools.utilities import smartLabels
import deeptools.getScorePerBigWigBin as score_bw
import deeptools.deepBlue as db

old_settings = np.seterr(all='ignore')


def parse_arguments(args=None):
    parser = \
        argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description="""

Given typically two or more bigWig files, ``multiBigwigSummary`` computes the average scores for each of the files in every genomic region.
This analysis is performed for the entire genome by running the program in ``bins`` mode, or for certain user selected regions in ``BED-file``
mode. Most commonly, the default output of ``multiBigwigSummary`` (a compressed numpy array, .npz) is used by other tools such as ``plotCorrelation`` or ``plotPCA`` for visualization and diagnostic purposes.

Note that using a single bigWig file is only recommended if you want to produce a bedGraph file (i.e., with the ``--outRawCounts`` option; the default output file cannot be used by ANY deepTools program if only a single file was supplied!).

A detailed sub-commands help is available by typing:

  multiBigwigSummary bins -h

  multiBigwigSummary BED-file -h


""",
            epilog='example usage:\n multiBigwigSummary bins '
                   '-b file1.bw file2.bw -o results.npz\n\n'
                   'multiBigwigSummary BED-file -b file1.bw file2.bw -o results.npz\n'
                   '--BED selection.bed'
                   ' \n\n',
            conflict_handler='resolve')

    parser.add_argument('--version', action='version',
                        version='multiBigwigSummary {}'.format(__version__))
    subparsers = parser.add_subparsers(
        title="commands",
        dest='command',
        metavar='')

    parent_parser = parserCommon.getParentArgParse(binSize=False)
    dbParser = parserCommon.deepBlueOptionalArgs()

    # bins mode options
    subparsers.add_parser(
        'bins',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[multiBigwigSummaryArgs(case='bins'),
                 parent_parser,
                 parserCommon.gtf_options(suppress=True),
                 dbParser
                 ],
        help="The average score is based on equally sized bins "
             "(10 kilobases by default), which consecutively cover the "
             "entire genome. The only exception is the last bin of a chromosome, which "
             "is often smaller. The output of this mode is commonly used to assess the "
             "overall similarity of different bigWig files.",
        add_help=False,
        usage='multiBigwigSummary bins '
              '-b file1.bw file2.bw '
              '-o results.npz\n')

    # BED file arguments
    subparsers.add_parser(
        'BED-file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[multiBigwigSummaryArgs(case='BED-file'),
                 parent_parser,
                 parserCommon.gtf_options(),
                 dbParser
                 ],
        help="The user provides a BED file that contains all regions "
             "that should be considered for the analysis. A "
             "common use is to compare scores (e.g. ChIP-seq scores) between "
             "different samples over a set of pre-defined peak regions.",
        usage='multiBigwigSummary BED-file '
              '-b file1.bw file2.bw '
              '-o results.npz --BED selection.bed\n',
        add_help=False)

    return parser


def process_args(args=None):
    args = parse_arguments().parse_args(args)

    if not args.labels and args.smartLabels:
        args.labels = smartLabels(args.bwfiles)
    elif not args.labels:
        args.labels = []
        for f in args.bwfiles:
            if f.startswith("http://") or f.startswith("https://") or f.startswith("ftp://"):
                args.labels.append(f.split("/")[-1])
            else:
                args.labels.append(os.path.basename(f))

    if len(args.bwfiles) != len(args.labels):
        sys.exit("The number of labels does not match the number of bigWig files.")

    return args


def multiBigwigSummaryArgs(case='bins'):
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    # define the arguments
    required.add_argument('--bwfiles', '-b',
                          metavar='FILE1 FILE2',
                          help='List of bigWig files, separated by spaces.',
                          nargs='+',
                          required=True)

    required.add_argument('--outFileName', '-out', '-o',
                          help='File name to save the compressed matrix file (npz format) '
                          'needed by the "plotPCA" and "plotCorrelation" tools.',
                          type=parserCommon.writableFile,
                          required=True)

    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")
    optional.add_argument('--labels', '-l',
                          metavar='sample1 sample2',
                          help='User defined labels instead of default labels from '
                          'file names. '
                          'Multiple labels have to be separated by spaces, e.g., '
                          '--labels sample1 sample2 sample3',
                          nargs='+')
    optional.add_argument('--smartLabels',
                          action='store_true',
                          help='Instead of manually specifying labels for the input '
                          'bigWig files, this causes deepTools to use the file name '
                          'after removing the path and extension.')

    optional.add_argument('--chromosomesToSkip',
                          metavar='chr1 chr2',
                          help='List of chromosomes that you do not want to be included. '
                          ' Useful to remove "random" or "extra" chr.',
                          nargs='+')

    if case == 'bins':
        optional.add_argument('--binSize', '-bs',
                              metavar='INT',
                              help='Size (in bases) of the windows sampled '
                              'from the genome. (Default: %(default)s)',
                              default=10000,
                              type=int)

        optional.add_argument('--distanceBetweenBins', '-n',
                              metavar='INT',
                              help='By default, multiBigwigSummary considers adjacent '
                              'bins of the specified --binSize. However, to '
                              'reduce the computation time, a larger distance '
                              'between bins can be given. Larger distances '
                              'results in fewer considered bins. (Default: %(default)s)',
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
                              help='Limits the analysis to '
                              'the regions specified in this file.',
                              metavar='file1.bed file2.bed',
                              nargs='+',
                              required=True)

    group = parser.add_argument_group('Output optional options')

    group.add_argument('--outRawCounts',
                       help='Save average scores per region for each bigWig file to a single tab-delimited file.',
                       type=parserCommon.writableFile,
                       metavar='FILE')

    return parser


def main(args=None):
    """
    1. get read counts at different positions either
    all of same length or from genomic regions from the BED file

    2. compute  the scores

    """
    args = process_args(args)

    if 'BED' in args:
        bed_regions = args.BED
    else:
        bed_regions = None

    if len(args.bwfiles) == 1 and not args.outRawCounts:
        sys.stderr.write("You've input a single bigWig file and not specified "
                         "--outRawCounts. The resulting output will NOT be "
                         "useful with any deepTools program!\n")

    # Preload deepBlue files, which need to then be deleted
    deepBlueFiles = []
    for idx, fname in enumerate(args.bwfiles):
        if db.isDeepBlue(fname):
            deepBlueFiles.append([fname, idx])
    if len(deepBlueFiles) > 0:
        sys.stderr.write("Preloading the following deepBlue files: {}\n".format(",".join([x[0] for x in deepBlueFiles])))
        if 'BED' in args:
            regs = db.makeRegions(args.BED, args)
        else:
            foo = db.deepBlue(deepBlueFiles[0][0], url=args.deepBlueURL, userKey=args.userKey)
            regs = db.makeTiles(foo, args)
            del foo
        for x in deepBlueFiles:
            x.extend([args, regs])
        if len(deepBlueFiles) > 1 and args.numberOfProcessors > 1:
            pool = multiprocessing.Pool(args.numberOfProcessors)
            res = pool.map_async(db.preloadWrapper, deepBlueFiles).get(9999999)
        else:
            res = list(map(db.preloadWrapper, deepBlueFiles))

        # substitute the file names with the temp files
        for (ftuple, r) in zip(deepBlueFiles, res):
            args.bwfiles[ftuple[1]] = r
        deepBlueFiles = [[x[0], x[1]] for x in deepBlueFiles]
        del regs

    num_reads_per_bin = score_bw.getScorePerBin(
        args.bwfiles,
        args.binSize,
        blackListFileName=args.blackListFileName,
        numberOfProcessors=args.numberOfProcessors,
        stepSize=args.binSize + args.distanceBetweenBins,
        verbose=args.verbose,
        region=args.region,
        bedFile=bed_regions,
        chrsToSkip=args.chromosomesToSkip,
        out_file_for_raw_data=args.outRawCounts,
        allArgs=args)

    sys.stderr.write("Number of bins "
                     "found: {}\n".format(num_reads_per_bin.shape[0]))

    if num_reads_per_bin.shape[0] < 2:
        exit("ERROR: too few non zero bins found.\n"
             "If using --region please check that this "
             "region is covered by reads.\n")

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
        f = open(args.outRawCounts, "r+")
        content = f.read()
        f.seek(0, 0)
        f.write(header + content)

        """
        if bed_regions:
            bed_regions.seek(0)
            reg_list = bed_regions.readlines()
            args.outRawCounts.write("#'chr'\t'start'\t'end'\t")
            args.outRawCounts.write("'" + "'\t'".join(args.labels) + "'\n")
            fmt = "\t".join(np.repeat('%s', num_reads_per_bin.shape[1])) + "\n"
            for idx, row in enumerate(num_reads_per_bin):
                args.outRawCounts.write("{}\t{}\t{}\t".format(*reg_list[idx].strip().split("\t")[0:3]))
                args.outRawCounts.write(fmt % tuple(row))

        else:
            args.outRawCounts.write("'" + "'\t'".join(args.labels) + "'\n")
            fmt = "\t".join(np.repeat('{}', num_reads_per_bin.shape[1])) + "\n"
            for row in num_reads_per_bin:
                args.outRawCounts.write(fmt.format(*tuple(row)))
        """
        f.close()

    # Clean up temporary bigWig files, if applicable
    if not args.deepBlueKeepTemp:
        for k, v in deepBlueFiles:
            os.remove(args.bwfiles[v])
    else:
        for k, v in deepBlueFiles:
            print("{} is stored in {}".format(k, args.bwfiles[v]))
