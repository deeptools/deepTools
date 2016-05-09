#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

# own tools
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

    parser.add_argument('--histogram', '-hist',
                        help='Save a .png file with a histogram '
                        'of the fragment length distribution.',
                        metavar='FILE')

    parser.add_argument('--numberOfProcessors', '-p',
                        help='Number of processors to use. The default is '
                        'to use 1.',
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
                        'the generated image. Leave blank for no title.',
                        default='')
    parser.add_argument('--maxFragmentLength',
                        help='The maximum fragment length in the histogram. A value of 0 (the default) indicates to use twice the mean fragment length',
                        default=0,
                        type=int)
    parser.add_argument('--logScale',
                        help='Plot on the log scale',
                        action='store_true')
    parser.add_argument('--binSize', '-bs',
                        metavar='INT',
                        help='Length in bases of the window used to sample the genome. (default 1000)',
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
                        'decreased. (default 1000000)',
                        default=1000000,
                        type=int)
    parser.add_argument('--blackListFileName', '-bl',
                        help="A BED file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered.",
                        metavar="BED file",
                        required=False)
    parser.add_argument('--verbose',
                        help='Set if processing data messages are wanted.',
                        action='store_true',
                        required=False)
    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def getFragSize(bam, args):
        fragment_len_dict, read_len_dict = get_read_and_fragment_length(bam, return_lengths=True,
                                                                        blackListFileName=args.blackListFileName,
                                                                        numberOfProcessors=args.numberOfProcessors,
                                                                        verbose=args.verbose,
                                                                        binSize=args.binSize,
                                                                        distanceBetweenBins=args.distanceBetweenBins)
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
        else:
            print("No pairs were found. Is the data from a paired-end sequencing experiment?")

        print("\nRead lengths:")
        print("Min.: {}\n1st Qu.: {}\nMean: {}\nMedian: {}\n"
              "3rd Qu.: {}\nMax.: {}\nStd: {}".format(read_len_dict['min'],
                                                      read_len_dict['qtile25'],
                                                      read_len_dict['mean'],
                                                      read_len_dict['median'],
                                                      read_len_dict['qtile75'],
                                                      read_len_dict['max'],
                                                      read_len_dict['std']))
        return fragment_len_dict


def main(args=None):
    args = parse_arguments().parse_args(args)

    fraglengths = {}
    for bam in args.bamfiles:
        fraglengths[bam] = getFragSize(bam, args)

    if args.histogram:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        if args.samplesLabel:
            if len(args.bamfiles) != len(args.samplesLabel):
                print("The number of labels does not match the number of BAM files.")
                exit(0)
            else:
                labels = args.samplesLabel
        else:
            labels = fraglengths.keys()

        i = 0
        for bam in fraglengths.keys():

            if args.maxFragmentLength > 0:
                maxVal = args.maxFragmentLength
            else:
                maxVal = fraglengths[bam]['mean'] * 2

            plt.hist(fraglengths[bam]['lengths'], 100,
                     range=(fraglengths[bam]['min'], maxVal),
                     alpha=0.5, label=labels[i],
                     log=args.logScale, normed=True)
            i += 1

        plt.xlabel('Fragment Length')
        plt.ylabel('Frequency')
        plt.legend(loc='upper right')
        plt.title(args.plotTitle)
        plt.savefig(args.histogram, bbox_inches=0)
        plt.close()


if __name__ == "__main__":
    main()
