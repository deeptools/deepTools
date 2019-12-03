#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
from deeptools import cm  # noqa: F401
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.stats import poisson

import plotly.offline as py
import plotly.graph_objs as go

import deeptools.countReadsPerBin as countR
import deeptools.sumCoveragePerBin as sumR
from deeptools import parserCommon
from deeptools.utilities import smartLabels

old_settings = np.seterr(all='ignore')
MAXLEN = 10000000


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

    if args.JSDsample is not None and args.JSDsample not in args.bamfiles:
        args.bamfiles.append(args.JSDsample)
        if args.labels and len(args.bamfiles) == len(args.labels) - 1:
            args.labels.append(args.JSDsample)

    if not args.labels:
        if args.smartLabels:
            args.labels = smartLabels(args.bamfiles)
        else:
            args.labels = args.bamfiles

    if len(args.bamfiles) != len(args.labels):
        sys.exit("The number of labels does not match the number of BAM files.")

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

    optional.add_argument('--smartLabels',
                          action='store_true',
                          help='Instead of manually specifying labels for the input '
                          'BAM/bigWig files, this causes deepTools to use the file name '
                          'after removing the path and extension.')

    optional.add_argument('--binSize', '-bs',
                          help='Window size in base pairs to '
                          'sample the genome. This times --numberOfSamples should be less than the genome size. (Default: %(default)s)',
                          default=500,
                          type=int)

    optional.add_argument('--numberOfSamples', '-n',
                          help='The number of bins that are sampled from the genome, '
                          'for which the overlapping number of reads is computed. (Default: %(default)s)',
                          default=5e5,
                          type=int)

    optional.add_argument('--plotFileFormat',
                          metavar='',
                          help='image format type. If given, this option '
                          'overrides the image format based on the ending '
                          'given via --plotFile '
                          'ending. The available options are: "png", '
                          '"eps", "pdf", "plotly" and "svg"',
                          choices=['png', 'pdf', 'svg', 'eps', 'plotly'])

    optional.add_argument('--plotTitle', '-T',
                          help='Title of the plot, to be printed on top of '
                          'the generated image. Leave blank for no title. (Default: %(default)s)',
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
                          'file and columns containing a number of metrics. '
                          'Please see the online documentation for a longer '
                          'explanation: http://deeptools.readthedocs.io/en/latest/content/feature/plotFingerprint_QC_metrics.html .',
                          type=parserCommon.writableFile,
                          metavar='FILE.txt')

    optional.add_argument('--JSDsample',
                          help='Reference sample against which to compute the '
                          'Jensen-Shannon distance and the CHANCE statistics. '
                          'If this is not specified, '
                          'then these will not be calculated. If '
                          '--outQualityMetrics is not specified then this will '
                          'be ignored. The Jensen-Shannon implementation is '
                          'based on code from Sitanshu Gakkhar at BCGSC. The '
                          'CHANCE implementation is based on code from Matthias '
                          'Haimel.',
                          metavar='sample.bam')

    return parser


def get_output_args():
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Output')
    group.add_argument('--plotFile', '-plot', '-o',
                       help='File name of the output figure. The file '
                       'ending will be used to determine the image '
                       'format. The available options are typically: "png", '
                       '"eps", "pdf" and "svg", e.g. : fingerprint.png.',
                       type=parserCommon.writableFile,
                       metavar='')

    group.add_argument('--outRawCounts',
                       help='Output file name to save the read counts per bin.',
                       type=parserCommon.writableFile,
                       metavar='')

    return parser


def binRelEntropy(p, q):
    """
    Return the relative binary entropy of x
    """
    x1 = 0
    x2 = 0
    if p > 0:
        x1 = p * np.log2(p / q)
    if p < 1:
        x2 = (1 - p) * np.log2((1 - p) / (1 - q))
    return np.fmax(0.0, x1 + x2)


def getCHANCE(args, idx, mat):
    """
    Compute the CHANCE p-value

    1) In short, sort IP from lowest to highest, cosorting input at the same time.
    2) Choose the argmax of the difference of the cumsum() of the above
    3) Determine a scale factor according to the ratio at the position at step 2.
    """
    # Get the index of the reference sample
    if args.JSDsample not in args.bamfiles:
        return ["NA", "NA", "NA"]
    refIdx = args.bamfiles.index(args.JSDsample)
    if refIdx == idx:
        return ["NA", "NA", "NA"]

    subMatrix = np.copy(mat[:, [idx, refIdx]])
    subMatrix[np.isnan(subMatrix)] = 0
    subMatrix = subMatrix[subMatrix[:, 0].argsort(), :]

    # Find the CHANCE statistic, which is the point of maximus difference
    cs = np.cumsum(subMatrix, axis=0)
    normed = cs / np.max(cs, axis=0).astype(float)
    csdiff = normed[:, 1] - normed[:, 0]
    k = np.argmax(csdiff)
    if csdiff[k] < 1e-6:
        # Don't bother with negative values
        return [0, 0, 0]
    p = normed[k, 0]  # Percent enrichment in IP
    q = normed[k, 1]  # Percent enrichment in input
    pcenrich = 100 * (len(csdiff) - k) / float(len(csdiff))
    diffenrich = 100.0 * (q - p)

    # CHANCE's JS divergence with binary entropy
    # Its p value is a ztest of this, which is largely useless IMO
    M = (p + q) / 2.0
    CHANCEdivergence = 0.5 * (binRelEntropy(p, M) + binRelEntropy(q, M))
    CHANCEdivergence = np.sqrt(CHANCEdivergence)

    return [pcenrich, diffenrich, CHANCEdivergence]


def getSyntheticJSD(vec):
    """
    This is largely similar to getJSD, with the 'input' sample being a Poisson distribution with lambda the average coverage in the IP bins
    """
    lamb = np.mean(vec)  # Average coverage
    coverage = np.sum(vec)

    chip = np.zeros(MAXLEN, dtype=np.int)
    for val in vec:
        # N.B., we need to clip past the end of the array
        if val >= MAXLEN:
            val = MAXLEN - 1
        # This effectively removes differences due to coverage percentages
        if val > 0:
            chip[int(val)] += 1
    input = coverage * poisson.pmf(np.arange(1, MAXLEN), lamb)
    if chip[-1] > 0:
        print("{} bins had coverage over the maximum value of {} during synthetic JSD computation".format(chip[-1], MAXLEN))

    return getJSDcommon(chip, input)


def getJSD(args, idx, mat):
    """
    Computes the Jensen-Shannon distance between two samples. This is essentially
    a symmetric version of Kullback-Leibler divergence. The implementation
    presented here is based on code from Sitanshu Gakkhar at BCGSC.

    Note that the interpolation has the effect of removing zero count coverage
    bins, which ends up being needed for the JSD calculation.

    args: The input arguments
    idx:  The column index of the current sample
    mat:  The matrix of counts
    """

    # Get the index of the reference sample
    if args.JSDsample not in args.bamfiles:
        return "NA"
    refIdx = args.bamfiles.index(args.JSDsample)
    if refIdx == idx:
        return "NA"

    # These will hold the coverage histograms
    chip = np.zeros(MAXLEN, dtype=np.int)
    input = np.zeros(MAXLEN, dtype=np.int)
    for row in mat:
        # ChIP
        val = row[idx]
        # N.B., we need to clip past the end of the array
        if val >= MAXLEN:
            val = MAXLEN - 1
        # This effectively removes differences due to coverage percentages
        if val > 0:
            chip[int(val)] += 1

        # Input
        val = row[refIdx]
        if val >= MAXLEN:
            val = MAXLEN - 1
        if val > 0:
            input[int(val)] += 1
    if input[-1] > 0:
        print("{} bins had coverage over the maximum value of {} in the input sample".format(input[-1], MAXLEN))
    if chip[-1] > 0:
        print("{} bins had coverage over the maximum value of {} in the ChIP sample".format(chip[-1], MAXLEN))

    return getJSDcommon(chip, input)


def getJSDcommon(chip, input):
    """
    This is a continuation of getJSD to allow getSyntheticJSD to reuse code
    """
    def signalAndBinDist(x):
        x = np.array(x)
        (n,) = x.shape
        signalValues = np.array(list(range(n)))
        totalSignal = x * signalValues
        normalizedTotalSignal = np.cumsum(totalSignal) / np.sum(totalSignal).astype("float")
        binDist = np.cumsum(x).astype("float") / sum(x)
        interpolater = interpolate.interp1d(binDist, normalizedTotalSignal, kind='linear', bounds_error=False, fill_value=(0, 1))
        return (binDist, normalizedTotalSignal, interpolater)

    # Interpolate the signals to evenly spaced bins, which also removes 0-coverage bins
    chipSignal = signalAndBinDist(chip)
    inputSignal = signalAndBinDist(input)

    # These are basically CDFs
    inputSignalInterp = inputSignal[2](np.arange(0, 1.00001, 0.00001))
    chipSignalInterp = chipSignal[2](np.arange(0, 1.00001, 0.00001))

    # If there are no low coverage bins then you can get nan as the first interpolated value.
    # That should instead be some small value
    if np.isnan(inputSignalInterp[0]):
        inputSignalInterp[0] = 1e-12
    if np.isnan(chipSignalInterp[0]):
        chipSignalInterp[0] = 1e-12

    # Differentiate to PMFs, do some sanity checking
    PMFinput = np.ediff1d(inputSignalInterp)
    PMFchip = np.ediff1d(chipSignalInterp)

    if abs(sum(PMFinput) - 1) > 0.01 or abs(sum(PMFchip) - 1) > 0.01:
        sys.stderr.write("Warning: At least one PMF integral is significantly different from 1! The JSD will not be returned")
        return "NA"

    # Compute the JSD from the PMFs
    M = (PMFinput + PMFchip) / 2.0
    JSD = 0.5 * (np.nansum(PMFinput * np.log2(PMFinput / M))) + 0.5 * (np.nansum(PMFchip * np.log2(PMFchip / M)))

    return np.sqrt(JSD)


def getExpected(mu):
    """
    Given a mean coverage mu, determine the AUC, X-intercept, and elbow point
    of a Poisson-distributed perfectly behaved input sample with the same coverage
    """
    x = np.arange(round(poisson.interval(0.99999, mu=mu)[1] + 1))  # This will be an appropriate range
    pmf = poisson.pmf(x, mu=mu)
    cdf = poisson.cdf(x, mu=mu)
    cs = np.cumsum(pmf * x)
    cs /= max(cs)
    XInt = cdf[np.nonzero(cs)[0][0]]
    AUC = sum(poisson.pmf(x, mu=mu) * cs)
    elbow = cdf[np.argmax(cdf - cs)]
    return (AUC, XInt, elbow)


def main(args=None):
    args = process_args(args)

    if not args.plotFile and not args.outRawCounts and not args.outQualityMetrics:
        sys.stderr.write("\nAt least one of --plotFile, --outRawCounts or --outQualityMetrics is required.\n")
        sys.exit(1)

    cr = sumR.SumCoveragePerBin(
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
        samFlag_exclude=args.samFlagExclude,
        minFragmentLength=args.minFragmentLength,
        maxFragmentLength=args.maxFragmentLength)

    num_reads_per_bin = cr.run()
    if num_reads_per_bin.sum() == 0:
        import sys
        sys.stderr.write(
            "\nNo reads were found in {} regions sampled. Check that the\n"
            "min mapping quality is not overly high and that the \n"
            "chromosome names between bam files are consistant.\n"
            "For small genomes, decrease the --numberOfSamples.\n"
            "\n".format(num_reads_per_bin.shape[0]))
        exit(1)

    if args.skipZeros:
        num_reads_per_bin = countR.remove_row_of_zeros(num_reads_per_bin)

    total = len(num_reads_per_bin[:, 0])
    x = np.arange(total).astype('float') / total  # normalize from 0 to 1

    if args.plotFile is not None:
        i = 0
        # matplotlib won't iterate through line styles by itself
        pyplot_line_styles = sum([7 * ["-"], 7 * ["--"], 7 * ["-."], 7 * [":"]], [])
        plotly_colors = ["#d73027", "#fc8d59", "#f33090", "#e0f3f8", "#91bfdb", "#4575b4"]
        plotly_line_styles = sum([6 * ["solid"], 6 * ["dot"], 6 * ["dash"], 6 * ["longdash"], 6 * ["dashdot"], 6 * ["longdashdot"]], [])
        data = []
        for i, reads in enumerate(num_reads_per_bin.T):
            count = np.cumsum(np.sort(reads))
            count = count / count[-1]  # to normalize y from 0 to 1
            if args.plotFileFormat == 'plotly':
                trace = go.Scatter(x=x, y=count, mode='lines', name=args.labels[i])
                trace['line'].update(dash=plotly_line_styles[i % 36], color=plotly_colors[i % 6])
                data.append(trace)
            else:
                j = i % len(pyplot_line_styles)
                plt.plot(x, count, label=args.labels[i], linestyle=pyplot_line_styles[j])
                plt.xlabel('rank')
                plt.ylabel('fraction w.r.t. bin with highest coverage')
        # set the plotFileFormat explicitly to None to trigger the
        # format from the file-extension
        if not args.plotFileFormat:
            args.plotFileFormat = None

        if args.plotFileFormat == 'plotly':
            fig = go.Figure()
            fig['data'] = data
            fig['layout'].update(title=args.plotTitle)
            fig['layout']['xaxis1'].update(title="rank")
            fig['layout']['yaxis1'].update(title="fraction w.r.t bin with highest coverage")
            py.plot(fig, filename=args.plotFile, auto_open=False)
        else:
            plt.legend(loc='upper left')
            plt.suptitle(args.plotTitle)
            plt.savefig(args.plotFile, bbox_inches=0, format=args.plotFileFormat)
            plt.close()

    if args.outRawCounts is not None:
        of = open(args.outRawCounts, "w")
        of.write("#plotFingerprint --outRawCounts\n")
        of.write("'" + "'\t'".join(args.labels) + "'\n")
        fmt = "\t".join(np.repeat('%d', num_reads_per_bin.shape[1])) + "\n"
        for row in num_reads_per_bin:
            of.write(fmt % tuple(row))
        of.close()

    if args.outQualityMetrics is not None:
        of = open(args.outQualityMetrics, "w")
        of.write("Sample\tAUC\tSynthetic AUC\tX-intercept\tSynthetic X-intercept\tElbow Point\tSynthetic Elbow Point")
        if args.JSDsample:
            of.write("\tJS Distance\tSynthetic JS Distance\t% genome enriched\tdiff. enrichment\tCHANCE divergence")
        else:
            of.write("\tSynthetic JS Distance")
        of.write("\n")
        line = np.arange(num_reads_per_bin.shape[0]) / float(num_reads_per_bin.shape[0] - 1)
        for idx, reads in enumerate(num_reads_per_bin.T):
            counts = np.cumsum(np.sort(reads))
            counts = counts / float(counts[-1])
            AUC = np.sum(counts) / float(len(counts))
            XInt = (np.argmax(counts > 0) + 1) / float(counts.shape[0])
            elbow = (np.argmax(line - counts) + 1) / float(counts.shape[0])
            expected = getExpected(np.mean(reads))  # A tuple of expected (AUC, XInt, elbow)
            of.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(args.labels[idx], AUC, expected[0], XInt, expected[1], elbow, expected[2]))
            if args.JSDsample:
                JSD = getJSD(args, idx, num_reads_per_bin)
                syntheticJSD = getSyntheticJSD(num_reads_per_bin[:, idx])
                CHANCE = getCHANCE(args, idx, num_reads_per_bin)
                of.write("\t{0}\t{1}\t{2}\t{3}\t{4}".format(JSD, syntheticJSD, CHANCE[0], CHANCE[1], CHANCE[2]))
            else:
                syntheticJSD = getSyntheticJSD(num_reads_per_bin[:, idx])
                of.write("\t{0}".format(syntheticJSD))
            of.write("\n")
        of.close()


if __name__ == "__main__":
    main()
