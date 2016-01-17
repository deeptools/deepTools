#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

# own tools
from deeptools.getFragmentAndReadSize import get_read_and_fragment_length
from deeptools._version import __version__


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Given a BAM file, several regions are sampled depending on the '
        'size of the genome and number of processors to estimate the paired-end '
        'fragment length. Properly paired reads are preferred for computation. ')
    parser.add_argument('bam',
                        help='BAM file to process',
                        metavar='bam-file')

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
    parser.add_argument('--plotTitle', '-T',
                        help='Title of the plot, to be printed on top of '
                        'the generated image. Leave blank for no title.',
                        default='')
    parser.add_argument('--verbose',
                        help='Set if processing data messages are wanted.',
                        action='store_true',
                        required=False)

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)
    fragment_len_dict, read_len_dict = get_read_and_fragment_length(args.bam, return_lengths=True,
                                                                    numberOfProcessors=args.numberOfProcessors,
                                                                    verbose=args.verbose)

    if fragment_len_dict['mean'] == 0:
        print "No pairs were found. Is the dataset paired-end?"

    print "Sample size: {}\n".format(fragment_len_dict['sample_size'])

    print "\nFragment lengths:"
    print "Min.: {}\n1st Qu.: {}\nMean: {}\nMedian: {}\n" \
          "3rd Qu.: {}\nMax.: {}\nStd: {}".format(fragment_len_dict['min'],
                                                  fragment_len_dict['qtile25'],
                                                  fragment_len_dict['mean'],
                                                  fragment_len_dict['median'],
                                                  fragment_len_dict['qtile75'],
                                                  fragment_len_dict['max'],
                                                  fragment_len_dict['std'])

    print "\nRead lengths:"
    print "Min.: {}\n1st Qu.: {}\nMean: {}\nMedian: {}\n" \
          "3rd Qu.: {}\nMax.: {}\nStd: {}".format(read_len_dict['min'],
                                                  read_len_dict['qtile25'],
                                                  read_len_dict['mean'],
                                                  read_len_dict['median'],
                                                  read_len_dict['qtile75'],
                                                  read_len_dict['max'],
                                                  read_len_dict['std'])

    if args.histogram:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.hist(fragment_len_dict['lengths'], 50,
                 range=(fragment_len_dict['min'], fragment_len_dict['mean'] * 2),
                 normed=True)
        plt.xlabel('Fragment Length')
        plt.ylabel('Frequency')
        plt.title(args.plotTitle)
        plt.savefig(args.histogram, bbox_inches=0)


if __name__ == "__main__":
    main()
