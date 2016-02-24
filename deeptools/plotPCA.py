#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from matplotlib import use as mplt_use
mplt_use('Agg')

from deeptools.correlation import Correlation
from deeptools._version import __version__


def parse_arguments(args=None):
    basic_args = plotCorrelationArgs()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Tool for generating a principal component analysis (PCA)
plot from multiBamSummary or multiBigwigSummary output.

Detailed help:

  plotPCA -h

""",
        epilog='example usages:\n'
               'plotPCA -in coverages.npz -o pca.png\n\n'
               ' \n\n',
        parents=[basic_args, ])
    return parser


def plotCorrelationArgs():
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    # define the arguments
    required.add_argument('--corData', '-in',
                          metavar='FILE',
                          help='Coverage file (generated by multiBamSummary or multiBigwigSummary)',
                          required=True)

    required.add_argument('--plotFile', '-o',
                          help='File name to save the plot to. '
                          'The extension determines the file format. '
                          'For example: '
                          'pca.pdf will save the PCA plot in PDF format. '
                          'The available options are: .png, '
                          '.eps, .pdf and .svg.',
                          type=argparse.FileType('w'),
                          metavar='FILE',
                          required=True)

    optional = parser.add_argument_group('Optional arguments')
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

    optional.add_argument('--plotFileFormat',
                          metavar='FILETYPE',
                          help='Image format type. If given, this option '
                          'overrides the image format based on the plotFile '
                          'ending. The available options are: png, '
                          'eps, pdf and svg.',
                          choices=['png', 'pdf', 'svg', 'eps'])

    optional.add_argument('--outFileNameData',
                          help='File name to save the data '
                          'underlying data for the average profile, e.g., '
                          'myProfile.tab.',
                          type=argparse.FileType('w'))

    optional.add_argument('--version', action='version',
                          version='%(prog)s {}'.format(__version__))

    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)

    corr = Correlation(args.corData,
                       labels=args.labels,)

    args.plotFile.close()

    corr.plot_pca(args.plotFile.name,
                  plot_title=args.plotTitle,
                  image_format=args.plotFileFormat)

    if args.outFileNameData is not None:
        import matplotlib
        mlab_pca = matplotlib.mlab.PCA(corr.matrix)
        n = len(corr.labels)
        of = args.outFileNameData
        of.write("Component\t{}\tEigenvalue\n".format("\t".join(corr.labels)))
        for i in xrange(n):
            of.write("{}".format(i + 1))
            for v in mlab_pca.Wt[i, :]:
                of.write("\t{}".format(v))
            of.write("\t{}\n".format(mlab_pca.s[i]))


if __name__ == "__main__":
    main()
