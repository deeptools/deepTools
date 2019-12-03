#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
from deeptools import cm  # noqa: F401
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import plotly.offline as py
import plotly.graph_objs as go

from deeptools.mapReduce import mapReduce, getUserRegion, blSubtract
from deeptools.getFragmentAndReadSize import get_read_and_fragment_length
from deeptools.utilities import getCommonChrNames, mungeChromosome, getTLen, smartLabels
from deeptools.bamHandler import openBam
from deeptoolsintervals import Enrichment, GTF
from deeptools.countReadsPerBin import CountReadsPerBin as cr
from deeptools import parserCommon


old_settings = np.seterr(all='ignore')


def parse_arguments(args=None):
    basic_args = plot_enrichment_args()

    # --region, --blackListFileName, -p and -v
    parent_parser = parserCommon.getParentArgParse(binSize=False)

    # --extend reads and such
    read_options = parserCommon.read_options()

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Tool for calculating and plotting the signal enrichment in either regions in BED
format or feature types (column 3) in GTF format. The underlying datapoints can also be output.
Metrics are plotted as a fraction of total reads. Regions in a BED file are assigned to the 'peak' feature.

detailed help:

  plotEnrichment -h

""",
        epilog='example usages:\n'
               'plotEnrichment -b file1.bam file2.bam --BED peaks.bed -o enrichment.png\n\n'
               ' \n\n',
        parents=[basic_args, parent_parser, read_options])

    return parser


def plot_enrichment_args():
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    # define the arguments
    required.add_argument('--bamfiles', '-b',
                          metavar='file1.bam file2.bam',
                          help='List of indexed bam files separated by spaces.',
                          nargs='+',
                          required=True)

    required.add_argument('--BED',
                          help='Limits the enrichment analysis to '
                          'the regions specified in these BED/GTF files. Enrichment '
                          'is calculated as the number of reads overlapping each '
                          'feature type. The feature type is column 3 in a GTF file '
                          'and "peak" for BED files.',
                          metavar='FILE1.bed FILE2.bed',
                          nargs='+',
                          required=True)

    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument('--plotFile', '-o',
                          help='File to save the plot to. The file extension determines the format, '
                          'so heatmap.pdf will save the heatmap in PDF format. '
                          'The available formats are: .png, '
                          '.eps, .pdf and .svg.',
                          type=parserCommon.writableFile,
                          metavar='FILE')

    optional.add_argument('--attributeKey',
                          help='Instead of deriving labels from the feature column in a GTF file, '
                          'use the given attribute key, such as gene_biotype. For BED files or '
                          'entries without the attribute key, None is used as the label.')

    optional.add_argument('--labels', '-l',
                          metavar='sample1 sample2',
                          help='User defined labels instead of default labels from '
                          'file names. '
                          'Multiple labels have to be separated by spaces, e.g. '
                          '--labels sample1 sample2 sample3',
                          nargs='+')

    optional.add_argument('--smartLabels',
                          action='store_true',
                          help='Instead of manually specifying labels for the input '
                          'BAM/BED/GTF files, this causes deepTools to use the file name '
                          'after removing the path and extension. For BED/GTF files, the '
                          'eventual region name will be overriden if specified inside '
                          'the file.')

    optional.add_argument('--regionLabels',
                          metavar="region1 region2",
                          help="For BED files, the label given to its region is "
                          "the file name, but this can be overridden by providing "
                          "a custom label. For GTF files this is ignored. Note "
                          "that if you provide labels, you MUST provide one for each "
                          "BED/GTF file, even though it will be ignored for GTF files.",
                          nargs='+')

    optional.add_argument('--plotTitle', '-T',
                          help='Title of the plot, to be printed on top of '
                          'the generated image. Leave blank for no title. (Default: %(default)s)',
                          default='')

    optional.add_argument('--plotFileFormat',
                          metavar='FILETYPE',
                          help='Image format type. If given, this option '
                          'overrides the image format based on the plotFile '
                          'ending. The available options are: png, '
                          'eps, pdf, plotly and svg.',
                          choices=['png', 'pdf', 'svg', 'eps', 'plotly'])

    optional.add_argument('--outRawCounts',
                          help='Save the counts per region to a tab-delimited file.',
                          type=parserCommon.writableFile,
                          metavar='FILE')

    optional.add_argument('--perSample',
                          help='Group the plots by sample, rather than by feature type (the default).',
                          action='store_true')

    optional.add_argument('--variableScales',
                          help='By default, the y-axes are always 0-100. This allows the axis range to be restricted.',
                          action='store_true')

    optional.add_argument('--plotHeight',
                          help='Plot height in cm. (Default: %(default)s)',
                          type=float,
                          default=20)

    optional.add_argument('--plotWidth',
                          help='Plot width in cm. The minimum value is 1 cm. (Default: %(default)s)',
                          type=float,
                          default=20)

    optional.add_argument('--colors',
                          help='List of colors to use '
                          'for the plotted lines. Color names '
                          'and html hex strings (e.g., #eeff22) '
                          'are accepted. The color names should '
                          'be space separated. For example, '
                          '--colors red blue green ',
                          nargs='+')

    optional.add_argument('--numPlotsPerRow',
                          help='Number of plots per row (Default: %(default)s)',
                          type=int,
                          default=4)

    optional.add_argument('--alpha',
                          default=0.9,
                          type=parserCommon.check_float_0_1,
                          help='The alpha channel (transparency) to use for the bars. '
                          'The default is 0.9 and values must be between 0 and 1.')

    optional.add_argument('--Offset',
                          help='Uses this offset inside of each read as the signal. This is useful in '
                          'cases like RiboSeq or GROseq, where the signal is 12, 15 or 0 bases past the '
                          'start of the read. This can be paired with the --filterRNAstrand option. '
                          'Note that negative values indicate offsets from the end of each read. A value '
                          'of 1 indicates the first base of the alignment (taking alignment orientation '
                          'into account). Likewise, a value of -1 is the last base of the alignment. An '
                          'offset of 0 is not permitted. If two values are specified, then they will be '
                          'used to specify a range of positions. Note that specifying something like '
                          '--Offset 5 -1 will result in the 5th through last position being used, which '
                          'is equivalent to trimming 4 bases from the 5-prime end of alignments.',
                          metavar='INT',
                          type=int,
                          nargs='+',
                          required=False)

    bed12 = parser.add_argument_group('BED12 arguments')

    bed12.add_argument('--keepExons',
                       help="For BED12 files, use each exon as a region, rather than columns 2/3",
                       action="store_true")

    return parser


def getBAMBlocks(read, defaultFragmentLength, centerRead, offset=None):
    """
    This is basically get_fragment_from_read from countReadsPerBin
    """
    blocks = None
    maxPairedFragmentLength = 0
    if defaultFragmentLength != "read length":
        maxPairedFragmentLength = 4 * defaultFragmentLength

    if defaultFragmentLength == 'read length':
        blocks = read.get_blocks()
    else:
        if cr.is_proper_pair(read, maxPairedFragmentLength):
            if read.is_reverse:
                fragmentStart = read.next_reference_start
                fragmentEnd = read.reference_end
            else:
                fragmentStart = read.reference_start
                # the end of the fragment is defined as
                # the start of the forward read plus the insert length
                fragmentEnd = read.reference_start + abs(read.template_length)
        # Extend using the default fragment length
        else:
            if read.is_reverse:
                fragmentStart = read.reference_end - defaultFragmentLength
                fragmentEnd = read.reference_end
            else:
                fragmentStart = read.reference_start
                fragmentEnd = read.reference_start + defaultFragmentLength
        if centerRead:
            fragmentCenter = fragmentEnd - (fragmentEnd - fragmentStart) / 2
            fragmentStart = fragmentCenter - read.infer_query_length(always=False) / 2
            fragmentEnd = fragmentStart + read.infer_query_length(always=False)

        assert fragmentStart < fragmentEnd, "fragment start greater than fragment" \
                                            "end for read {}".format(read.query_name)
        blocks = [(int(fragmentStart), int(fragmentEnd))]

    # Handle read offsets, if needed
    if offset is not None:
        rv = [(None, None)]
        if len(offset) > 1:
            if offset[0] > 0:
                offset[0] -= 1
            if offset[1] < 0:
                offset[1] += 1
        else:
            if offset[0] > 0:
                offset[0] -= 1
                offset = [offset[0], offset[0] + 1]
            else:
                offset = [offset[0], None]
        if offset[1] == 0:
            # -1 gets switched to 0, which screws things up
            offset = (offset[0], None)

        stretch = []
        # For the sake of simplicity, convert [(10, 20), (30, 40)] to [10, 11, 12, 13, ..., 40]
        # Then subset accordingly
        for block in blocks:
            stretch.extend(range(block[0], block[1]))
        if read.is_reverse:
            stretch = stretch[::-1]
        try:
            foo = stretch[offset[0]:offset[1]]
        except:
            return rv

        if len(foo) == 0:
            return rv
        if read.is_reverse:
            foo = foo[::-1]
        # Convert the stretch back to a list of tuples
        foo = np.array(foo)
        d = foo[1:] - foo[:-1]
        idx = np.argwhere(d > 1).flatten().tolist()  # This now holds the interval bounds as a list
        idx.append(-1)
        last = 0
        blocks = []
        for i in idx:
            blocks.append((foo[last].astype("int"), foo[i].astype("int") + 1))
            last = i + 1
    return blocks


def getEnrichment_worker(arglist):
    """
    This is the worker function of plotEnrichment.

    In short, given a region, iterate over all reads **starting** in it.
    Filter/extend them as requested and check each for an overlap with
    findOverlaps. For each overlap, increment the counter for that feature.
    """
    chrom, start, end, args, defaultFragmentLength = arglist
    if args.verbose:
        sys.stderr.write("Processing {}:{}-{}\n".format(chrom, start, end))

    olist = []
    total = [0] * len(args.bamfiles)
    for idx, f in enumerate(args.bamfiles):
        odict = dict()
        for x in gtf.features:
            odict[x] = 0
        fh = openBam(f)

        chrom = mungeChromosome(chrom, fh.references)

        lpos = None
        prev_pos = set()
        for read in fh.fetch(chrom, start, end):
            # Filter
            if read.pos < start:
                # Ensure that a given alignment is processed only once
                continue
            if read.flag & 4:
                continue
            if args.minMappingQuality and read.mapq < args.minMappingQuality:
                continue
            if args.samFlagInclude and read.flag & args.samFlagInclude != args.samFlagInclude:
                continue
            if args.samFlagExclude and read.flag & args.samFlagExclude != 0:
                continue
            tLen = getTLen(read)
            if args.minFragmentLength > 0 and tLen < args.minFragmentLength:
                continue
            if args.maxFragmentLength > 0 and tLen > args.maxFragmentLength:
                continue
            if args.ignoreDuplicates:
                # Assuming more or less concordant reads, use the fragment bounds, otherwise the start positions
                if tLen >= 0:
                    s = read.pos
                    e = s + tLen
                else:
                    s = read.pnext
                    e = s - tLen
                if read.reference_id != read.next_reference_id:
                    e = read.pnext
                if lpos is not None and lpos == read.reference_start \
                        and (s, e, read.next_reference_id, read.is_reverse) in prev_pos:
                    continue
                if lpos != read.reference_start:
                    prev_pos.clear()
                lpos = read.reference_start
                prev_pos.add((s, e, read.next_reference_id, read.is_reverse))
            total[idx] += 1

            # Get blocks, possibly extending
            features = gtf.findOverlaps(chrom, getBAMBlocks(read, defaultFragmentLength, args.centerReads, args.Offset))

            if features is not None and len(features) > 0:
                for x in features:
                    odict[x] += 1
        olist.append(odict)
    return olist, gtf.features, total


def plotEnrichment(args, featureCounts, totalCounts, features):
    # get the number of rows and columns
    if args.perSample:
        totalPlots = len(args.bamfiles)
        barsPerPlot = len(features)
    else:
        totalPlots = len(features)
        barsPerPlot = len(args.bamfiles)
    cols = min(args.numPlotsPerRow, totalPlots)
    rows = np.ceil(totalPlots / float(args.numPlotsPerRow)).astype(int)

    # Handle the colors
    if not args.colors:
        cmap_plot = plt.get_cmap('jet')
        args.colors = cmap_plot(np.arange(barsPerPlot, dtype=float) / float(barsPerPlot))
        if args.plotFileFormat == 'plotly':
            args.colors = range(barsPerPlot)
    elif len(args.colors) < barsPerPlot:
        sys.exit("Error: {0} colors were requested, but {1} were needed!".format(len(args.colors), barsPerPlot))

    data = []
    if args.plotFileFormat == 'plotly':
        fig = go.Figure()
        fig['layout'].update(title=args.plotTitle)
        domainWidth = .9 / cols
        domainHeight = .9 / rows
        bufferHeight = 0.0
        if rows > 1:
            bufferHeight = 0.1 / (rows - 1)
        bufferWidth = 0.0
        if cols > 1:
            bufferWidth = 0.1 / (cols - 1)
    else:
        grids = gridspec.GridSpec(rows, cols)
        plt.rcParams['font.size'] = 10.0

        # convert cm values to inches
        fig = plt.figure(figsize=(args.plotWidth / 2.54, args.plotHeight / 2.54))
        fig.suptitle(args.plotTitle, y=(1 - (0.06 / args.plotHeight)))

    for i in range(totalPlots):
        col = i % cols
        row = np.floor(i / float(args.numPlotsPerRow)).astype(int)

        if args.perSample:
            xlabels = features
            ylabel = "% alignments in {0}".format(args.labels[i])
            vals = [featureCounts[i][foo] for foo in features]
            vals = 100 * np.array(vals, dtype='float64') / totalCounts[i]
        else:
            xlabels = args.labels
            ylabel = "% {0}".format(features[i])
            vals = [foo[features[i]] for foo in featureCounts]
            vals = 100 * np.array(vals, dtype='float64') / np.array(totalCounts, dtype='float64')

        if args.plotFileFormat == 'plotly':
            xanchor = 'x{}'.format(i + 1)
            yanchor = 'y{}'.format(i + 1)
            base = row * (domainHeight + bufferHeight)
            domain = [base, base + domainHeight]
            fig['layout']['xaxis{}'.format(i + 1)] = {'domain': domain, 'anchor': yanchor}
            base = col * (domainWidth + bufferWidth)
            domain = [base, base + domainWidth]
            fig['layout']['yaxis{}'.format(i + 1)] = {'domain': domain, 'anchor': xanchor, 'title': ylabel}
            if args.variableScales is False:
                fig['layout']['yaxis{}'.format(i + 1)].update(range=[0, 100])
            trace = go.Bar(x=xlabels,
                           y=vals,
                           opacity=args.alpha,
                           orientation='v',
                           showlegend=False,
                           xaxis=xanchor,
                           yaxis=yanchor,
                           name=ylabel,
                           marker={'color': args.colors, 'line': {'color': args.colors}})
            data.append(trace)
        else:
            ax = plt.subplot(grids[row, col])
            ax.bar(np.arange(vals.shape[0]), vals, width=1.0, bottom=0.0, align='center', color=args.colors, edgecolor=args.colors, alpha=args.alpha)
            ax.set_ylabel(ylabel)
            ax.set_xticks(np.arange(vals.shape[0]))
            ax.set_xticklabels(xlabels, rotation='vertical')
            if args.variableScales is False:
                ax.set_ylim(0.0, 100.0)

    if args.plotFileFormat == 'plotly':
        fig['data'] = data
        py.plot(fig, filename=args.plotFile, auto_open=False)
        # colors
    else:
        plt.subplots_adjust(wspace=0.05, hspace=0.3, bottom=0.15, top=0.80)
        plt.tight_layout()
        plt.savefig(args.plotFile, dpi=200, format=args.plotFileFormat)
        plt.close()


def getChunkLength(args, chromSize):
    """
    There's no point in parsing the GTF time over and over again needlessly.
    Emprically, it seems that adding ~4x the number of workers is ideal, since
    coverage is non-uniform. This is a heuristic way of approximating that.

    Note that if there are MANY small contigs and a few large ones (e.g., the
    max and median lengths are >10x different, then it's best to take a
    different tack.
    """

    if args.region:
        chromSize, region_start, region_end, genomeChunkLength = getUserRegion(chromSize, args.region)
        rv = np.ceil((region_start - region_end) / float(4 * args.numberOfProcessors)).astype(int)
        return max(1, rv)

    bl = None
    if args.blackListFileName:
        bl = GTF(args.blackListFileName)

    lengths = []
    for k, v in chromSize:
        regs = blSubtract(bl, k, [0, v])
        for reg in regs:
            lengths.append(reg[1] - reg[0])

    if len(lengths) >= 4 * args.numberOfProcessors:
        rv = np.median(lengths).astype(int)
        # In cases like dm6 or GRCh38, there are a LOT of really small contigs, which will cause the median to be small and performance to tank
        if np.max(lengths) >= 10 * rv:
            rv = np.ceil(np.sum(lengths) / (4.0 * args.numberOfProcessors)).astype(int)
    else:
        rv = np.ceil(np.sum(lengths) / (4.0 * args.numberOfProcessors)).astype(int)

    return max(1, rv)


def main(args=None):

    args = parse_arguments().parse_args(args)

    if not args.outRawCounts and not args.plotFile:
        sys.exit("Error: You need to specify at least one of --plotFile or --outRawCounts!\n")

    if args.labels is None:
        args.labels = args.bamfiles
    if args.smartLabels:
        args.labels = smartLabels(args.bamfiles)
    if len(args.labels) != len(args.bamfiles):
        sys.exit("Error: The number of labels ({0}) does not match the number of BAM files ({1})!".format(len(args.labels), len(args.bamfiles)))

    # Ensure that if we're given an attributeKey that it's not empty
    if args.attributeKey and args.attributeKey == "":
        args.attributeKey = None

    global gtf
    if not args.regionLabels and args.smartLabels:
        args.regionLabels = smartLabels(args.BED)
    gtf = Enrichment(args.BED, keepExons=args.keepExons, labels=args.regionLabels, attributeKey=args.attributeKey)

    # Get fragment size and chromosome dict
    fhs = [openBam(x) for x in args.bamfiles]
    chromSize, non_common_chr = getCommonChrNames(fhs, verbose=args.verbose)
    for fh in fhs:
        fh.close()

    frag_len_dict, read_len_dict = get_read_and_fragment_length(args.bamfiles[0],
                                                                return_lengths=False,
                                                                blackListFileName=args.blackListFileName,
                                                                numberOfProcessors=args.numberOfProcessors,
                                                                verbose=args.verbose)
    if args.extendReads:
        if args.extendReads is True:
            # try to guess fragment length if the bam file contains paired end reads
            if frag_len_dict:
                defaultFragmentLength = frag_len_dict['median']
            else:
                sys.exit("*ERROR*: library is not paired-end. Please provide an extension length.")
            if args.verbose:
                print("Fragment length based on paired en data "
                      "estimated to be {0}".format(frag_len_dict['median']))
        elif args.extendReads < read_len_dict['median']:
            sys.stderr.write("*WARNING*: read extension is smaller than read length (read length = {}). "
                             "Reads will not be extended.\n".format(int(read_len_dict['median'])))
            defaultFragmentLength = 'read length'
        elif args.extendReads > 2000:
            sys.exit("*ERROR*: read extension must be smaller that 2000. Value give: {} ".format(args.extendReads))
        else:
            defaultFragmentLength = args.extendReads
    else:
        defaultFragmentLength = 'read length'

    # Get the chunkLength
    chunkLength = getChunkLength(args, chromSize)

    # Map reduce to get the counts/file/feature
    res = mapReduce([args, defaultFragmentLength],
                    getEnrichment_worker,
                    chromSize,
                    genomeChunkLength=chunkLength,
                    region=args.region,
                    blackListFileName=args.blackListFileName,
                    numberOfProcessors=args.numberOfProcessors,
                    verbose=args.verbose)

    features = res[0][1]
    featureCounts = []
    for i in list(range(len(args.bamfiles))):
        d = dict()
        for x in features:
            d[x] = 0
        featureCounts.append(d)

    # res is a list, with each element a list (length len(args.bamfiles)) of dicts
    totalCounts = [0] * len(args.bamfiles)
    for x in res:
        for i, y in enumerate(x[2]):
            totalCounts[i] += y
        for i, y in enumerate(x[0]):
            for k, v in y.items():
                featureCounts[i][k] += v

    # Make a plot
    if args.plotFile:
        plotEnrichment(args, featureCounts, totalCounts, features)

    # Raw counts
    if args.outRawCounts:
        of = open(args.outRawCounts, "w")
        of.write("file\tfeatureType\tpercent\tfeatureReadCount\ttotalReadCount\n")
        for i, x in enumerate(args.labels):
            for k, v in featureCounts[i].items():
                of.write("{0}\t{1}\t{2:5.2f}\t{3}\t{4}\n".format(x, k, (100.0 * v) / totalCounts[i], v, totalCounts[i]))
        of.close()
