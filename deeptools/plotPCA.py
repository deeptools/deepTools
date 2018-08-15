#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

from deeptools.correlation import Correlation
from deeptools._version import __version__


def parse_arguments(args=None):
    basic_args = plotCorrelationArgs()
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Tool for generating a principal component analysis (PCA)
plot from multiBamSummary or multiBigwigSummary output. By default, the loadings for each sample in each principal component is plotted. If the data is transposed, the projections of each sample on the requested principal components is plotted instead.

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

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--plotFile', '-o',
                          help='File name to save the plot to. '
                          'The extension determines the file format. '
                          'For example: '
                          'pca.pdf will save the PCA plot in PDF format. '
                          'The available options are: .png, '
                          '.eps, .pdf and .svg. If this option is omitted, then you MUST specify --outFileNameData',
                          metavar='FILE')

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
                          'eps, pdf, plotly and svg.',
                          choices=['png', 'pdf', 'svg', 'eps', 'plotly'])

    optional.add_argument('--plotHeight',
                          help='Plot height in cm.',
                          type=float,
                          default=10)

    optional.add_argument('--plotWidth',
                          help='Plot width in cm. The minimum value is 1 cm.',
                          type=float,
                          default=10)

    optional.add_argument('--outFileNameData',
                          metavar='file.tab',
                          help='File name to which the data underlying the plot '
                          'should be saved, such as myPCA.tab. For untransposed '
                          'data, this is the loading per-sample and PC as well '
                          'as the eigenvalues. For transposed data, this is the '
                          'rotation per-sample and PC and the eigenvalues. The '
                          'projections are truncated to the number of '
                          'eigenvalues for transposed data.')

    optional.add_argument('--ntop',
                          help='Use only the top N most variable rows in the '
                          'original matrix. Specifying 0 will result in all '
                          'rows being used. If the matrix is to be transposed, '
                          'rows with 0 variance are always excluded, even if a '
                          'values of 0 is specified. The default is 1000.',
                          type=int,
                          default=1000)

    optional.add_argument('--PCs',
                          help='The principal components to plot. If specified, '
                          'you must provide two different integers, greater '
                          'than zero, separated by a space. An example (and the default) is "1 2".',
                          type=int,
                          nargs=2,
                          default=[1, 2])

    optional.add_argument('--log2',
                          help='log2 transform the datapoints prior to computing '
                          'the PCA. Note that 0.01 is added to all values to '
                          'prevent 0 values from becoming -infinity. Using this '
                          'option with input that contains negative values will '
                          'result in an error.',
                          action='store_true')

    optional.add_argument('--colors',
                          metavar="COLORS",
                          nargs='+',
                          help="A list of colors for the symbols. Color names and html hex string (e.g., #eeff22) are accepted. The color names should be space separated. For example, --colors red blue green. If not specified, the symbols will be given automatic colors.")

    optional.add_argument('--markers',
                          metavar="MARKERS",
                          nargs='+',
                          help="A list of markers for the symbols. (e.g., '<','>','o') are accepted. The marker values should be space separated. For example, --markers 's' 'o' 's' 'o'. If not specified, the symbols will be given automatic shapes.")

    optional.add_argument('--version', action='version',
                          version='%(prog)s {}'.format(__version__))

    optionalEx = optional.add_mutually_exclusive_group()
    optionalEx.add_argument('--transpose',
                            help='Perform the PCA on the transposed matrix, (i.e., on the '
                            'matrix where rows are samples and columns are '
                            'bins/features. This then matches what is typically '
                            'done in R.',
                            action='store_true')

    optionalEx.add_argument('--rowCenter',
                            help='When specified, each row (bin, gene, etc.) '
                            'in the matrix is centered at 0 before the PCA is '
                            'computed. This is useful only if you have a strong '
                            'bin/gene/etc. correlation and the resulting '
                            'principal component has samples stacked vertically. This option is not applicable if --transpose is specified.',
                            action='store_true')

    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)

    if args.plotFile is None and args.outFileNameData is None:
        sys.exit("At least one of --plotFile and --outFileNameData must be specified!\n")

    if args.ntop < 0:
        sys.exit("The value specified for --ntop must be >= 0!\n")

    if args.PCs[0] == args.PCs[1]:
        sys.exit("You must specify different principal components!\n")
    if args.PCs[0] <= 0 or args.PCs[1] <= 0:
        sys.exit("The specified principal components must be at least 1!\n")

    corr = Correlation(args.corData,
                       labels=args.labels,)

    corr.rowCenter = args.rowCenter
    corr.transpose = args.transpose
    corr.ntop = args.ntop
    corr.log2 = args.log2

    Wt, eigenvalues = corr.plot_pca(args.plotFile,
                                    PCs=args.PCs,
                                    plot_title=args.plotTitle,
                                    image_format=args.plotFileFormat,
                                    plotWidth=args.plotWidth,
                                    plotHeight=args.plotHeight,
                                    cols=args.colors,
                                    marks=args.markers)

    if args.outFileNameData is not None:
        of = open(args.outFileNameData, "w")
        of.write("#plotPCA --outFileNameData\n")
        of.write("Component\t{}\tEigenvalue\n".format("\t".join(corr.labels)))
        n = eigenvalues.shape[0]
        for i in range(n):
            of.write("{}\t{}\t{}\n".format(i + 1, "\t".join(["{}".format(x) for x in Wt[i, :]]), eigenvalues[i]))
        of.close()


if __name__ == "__main__":
    main()
