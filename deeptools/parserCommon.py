import argparse
import deeptools.config as cfg
import os
from deeptools._version import __version__


def check_float_0_1(value):
    v = float(value)
    if v < 0.0 or v > 1.0:
        raise argparse.ArgumentTypeError("%s is an invalid floating point value. It must be between 0.0 and 1.0" % value)
    return v


def output(args=None):
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Output')
    group.add_argument('--outFileName', '-o',
                       help='Output file name.',
                       metavar='FILENAME',
                       type=writableFile,
                       required=True)

    group.add_argument('--outFileFormat', '-of',
                       help='Output file type. Either "bigwig" or "bedgraph".',
                       choices=['bigwig', 'bedgraph'],
                       default='bigwig')

    return parser


def read_options():
    """Common arguments related to BAM files and the interpretation
    of the read coverage
    """
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Read processing options')

    group.add_argument('--extendReads', '-e',
                       help='This parameter allows the extension of reads to '
                       'fragment size. If set, each read is '
                       'extended, without exception.\n'
                       '*NOTE*: This feature is generally NOT recommended for '
                       'spliced-read data, such as RNA-seq, as it would '
                       'extend reads over skipped regions.\n'
                       '*Single-end*: Requires a user specified value for the '
                       'final fragment length. Reads that already exceed this '
                       'fragment length will not be extended.\n'
                       '*Paired-end*: Reads with mates are always extended to '
                       'match the fragment size defined by the two read mates. '
                       'Unmated reads, mate reads that map too far apart '
                       '(>4x fragment length) or even map to different '
                       'chromosomes are treated like single-end reads. The input '
                       'of a fragment length value is optional. If '
                       'no value is specified, it is estimated from the '
                       'data (mean of the fragment size of all mate reads).\n',
                       type=int,
                       nargs='?',
                       const=True,
                       default=False,
                       metavar="INT bp")

    group.add_argument('--ignoreDuplicates',
                       help='If set, reads that have the same orientation '
                       'and start position will be considered only '
                       'once. If reads are paired, the mate\'s position '
                       'also has to coincide to ignore a read.',
                       action='store_true'
                       )

    group.add_argument('--minMappingQuality',
                       metavar='INT',
                       help='If set, only reads that have a mapping '
                       'quality score of at least this are '
                       'considered.',
                       type=int,
                       )

    group.add_argument('--centerReads',
                       help='By adding this option, reads are centered with '
                       'respect to the fragment length. For paired-end data, '
                       'the read is centered at the fragment length defined '
                       'by the two ends of the fragment. For single-end data, the '
                       'given fragment length is used. This option is '
                       'useful to get a sharper signal around enriched '
                       'regions.',
                       action='store_true')

    group.add_argument('--samFlagInclude',
                       help='Include reads based on the SAM flag. For example, '
                       'to get only reads that are the first mate, use a flag of 64. '
                       'This is useful to count properly paired reads only once, '
                       'as otherwise the second mate will be also considered for the '
                       'coverage.',
                       metavar='INT',
                       default=None,
                       type=int,
                       required=False)

    group.add_argument('--samFlagExclude',
                       help='Exclude reads based on the SAM flag. For example, '
                       'to get only reads that map to the forward strand, use '
                       '--samFlagExclude 16, where 16 is the SAM flag for reads '
                       'that map to the reverse strand.',
                       metavar='INT',
                       default=None,
                       type=int,
                       required=False)

    group.add_argument('--minFragmentLength',
                       help='The minimum fragment length needed for read/pair '
                       'inclusion. Note that a value other than 0 will exclude '
                       'all single-end reads. This option is primarily useful '
                       'in ATACseq experiments, for filtering mono- or '
                       'di-nucleosome fragments.',
                       metavar='INT',
                       default=0,
                       type=int,
                       required=False)

    group.add_argument('--maxFragmentLength',
                       help='The maximum fragment length needed for read/pair '
                       'inclusion. A value of 0 disables filtering and is '
                       'needed for including single-end and orphan reads.',
                       metavar='INT',
                       default=0,
                       type=int,
                       required=False)

    return parser


def gtf_options(suppress=False):
    """
    Arguments present whenever a BED/GTF file can be used
    """
    if suppress:
        parser = argparse.ArgumentParser(add_help=False)
        group = parser
    else:
        parser = argparse.ArgumentParser(add_help=False)
        group = parser.add_argument_group('GTF/BED12 options')

    if suppress:
        help = argparse.SUPPRESS
    else:
        help = 'When either a BED12 or GTF file are used to provide \
        regions, perform the computation on the merged exons, \
        rather than using the genomic interval defined by the \
        5-prime and 3-prime most transcript bound (i.e., columns \
        2 and 3 of a BED file). If a BED3 or BED6 file is used \
        as input, then columns 2 and 3 are used as an exon.'

    group.add_argument('--metagene',
                       help=help,
                       action='store_true',
                       dest='keepExons')

    if suppress is False:
        help = 'When a GTF file is used to provide regions, only \
        entries with this value as their feature (column 2) \
        will be processed as transcripts.'

    group.add_argument('--transcriptID',
                       help=help,
                       default='transcript')

    if suppress is False:
        help = 'When a GTF file is used to provide regions, only \
        entries with this value as their feature (column 2) \
        will be processed as exons. CDS would be another common \
        value for this.'

    group.add_argument('--exonID',
                       help=help,
                       default='exon')

    if suppress is False:
        help = 'Each region has an ID (e.g., ACTB) assigned to it, \
        which for BED files is either column 4 (if it exists) \
        or the interval bounds. For GTF files this is instead \
        stored in the last column as a key:value pair (e.g., as \
        \'transcript_id "ACTB"\', for a key of transcript_id \
        and a value of ACTB). In some cases it can be \
        convenient to use a different identifier. To do so, set \
        this to the desired key.'

    group.add_argument('--transcript_id_designator',
                       help=help,
                       default='transcript_id')

    return parser


def normalization_options():
    """Common arguments related to read coverage normalization
    """
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Read coverage normalization options')

    group.add_argument('--normalizeTo1x',
                       help='Report read coverage normalized to 1x '
                       'sequencing depth (also known as Reads Per Genomic '
                       'Content (RPGC)). Sequencing depth is defined as: '
                       '(total number of mapped reads * fragment length) / '
                       'effective genome size.\nThe scaling factor used '
                       'is the inverse of the sequencing depth computed '
                       'for the sample to match the 1x coverage. '
                       'To use this option, the '
                       'effective genome size has to be indicated after the '
                       'option. The effective genome size is the portion '
                       'of the genome that is mappable. Large fractions of '
                       'the genome are stretches of NNNN that should be '
                       'discarded. Also, if repetitive regions were not '
                       'included in the mapping of reads, the effective '
                       'genome size needs to be adjusted accordingly. '
                       'Common values are: mm9: 2,150,570,000; '
                       'hg19:2,451,960,000; dm3:121,400,000 and ce10:93,260,000. '
                       'See Table 2 of http://www.plosone.org/article/info:doi/10.1371/journal.pone.0030377 '
                       'or http://www.nature.com/nbt/journal/v27/n1/fig_tab/nbt.1518_T1.html '
                       'for several effective genome sizes.',
                       metavar='EFFECTIVE GENOME SIZE LENGTH',
                       default=None,
                       type=int,
                       required=False)

    group.add_argument('--normalizeUsingRPKM',
                       help='Use Reads Per Kilobase per Million reads to '
                       'normalize the number of reads per bin. The formula '
                       'is: RPKM (per bin) =  number of reads per bin / '
                       '( number of mapped reads (in millions) * bin '
                       'length (kb) ). Each read is considered independently,'
                       'if you want to only count either of the mate pairs in'
                       'paired-end data, use the --samFlag option.',
                       action='store_true',
                       required=False)

    group.add_argument('--ignoreForNormalization', '-ignore',
                       help='A list of space-delimited chromosome names '
                       'containing those chromosomes that should be excluded '
                       'for computing the normalization. This is useful when considering '
                       'samples with unequal coverage across chromosomes, like male '
                       'samples. An usage examples is --ignoreForNormalization chrX chrM.',
                       nargs='+')

    group.add_argument('--skipNonCoveredRegions', '--skipNAs',
                       help='This parameter determines if non-covered regions '
                       '(regions without overlapping reads) in a BAM file should '
                       'be skipped. The default is to treat those regions as having a value of zero. '
                       'The decision to skip non-covered regions '
                       'depends on the interpretation of the data. Non-covered regions '
                       'may represent, for example, repetitive regions that should be skipped.',
                       action='store_true')

    group.add_argument('--smoothLength',
                       metavar="INT bp",
                       help='The smooth length defines a window, larger than '
                       'the binSize, to average the number of reads. For '
                       'example, if the --binSize is set to 20 and the '
                       '--smoothLength is set to 60, then, for each '
                       'bin, the average of the bin and its left and right '
                       'neighbors is considered. Any value smaller than '
                       '--binSize will be ignored and no smoothing will be '
                       'applied.',
                       type=int)

    return parser


def getParentArgParse(args=None, binSize=True, blackList=True):
    """
    Typical arguments for several tools
    """

    parser = argparse.ArgumentParser(add_help=False)
    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument('--version', action='version',
                          version='%(prog)s {}'.format(__version__))

    if binSize:
        optional.add_argument('--binSize', '-bs',
                              help='Size of the bins, in bases, for the output '
                              'of the bigwig/bedgraph file.',
                              metavar="INT bp",
                              type=int,
                              default=50)

    optional.add_argument('--region', '-r',
                          help='Region of the genome to limit the operation '
                          'to - this is useful when testing parameters to '
                          'reduce the computing time. The format is '
                          'chr:start:end, for example --region chr10 or '
                          '--region chr10:456700:891000.',
                          metavar="CHR:START:END",
                          required=False,
                          type=genomicRegion)

    if blackList:
        optional.add_argument('--blackListFileName', '-bl',
                              help="A BED or GTF file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant.",
                              metavar="BED file",
                              nargs="+",
                              required=False)

    optional.add_argument('--numberOfProcessors', '-p',
                          help='Number of processors to use. Type "max/2" to '
                          'use half the maximum number of processors or "max" '
                          'to use all available processors.',
                          metavar="INT",
                          type=numberOfProcessors,
                          default=cfg.config.get('general',
                                                 'default_proc_number'),
                          required=False)

    optional.add_argument('--verbose', '-v',
                          help='Set to see processing messages.',
                          action='store_true')

    return parser


def numberOfProcessors(string):
    import multiprocessing
    availProc = multiprocessing.cpu_count()

    if string == "max/2":  # default case
        # by default half of the available processors are used
        numberOfProcessors = int(availProc * 0.5)
    elif string == "max":
        # use all available processors
        numberOfProcessors = availProc
    else:
        try:
            numberOfProcessors = int(string)
        except ValueError:
            raise argparse.ArgumentTypeError(
                "{} is not a valid number of processors".format(string))

        except Exception as e:
            raise argparse.ArgumentTypeError("the value given is not valid. "
                                             "Error message: {}\nThe number of "
                                             "available processors in your "
                                             "computer is {}.".format(string, e,
                                                                      availProc))

        if numberOfProcessors > availProc:
            numberOfProcessors = availProc

    return numberOfProcessors


def genomicRegion(string):
    # remove whitespaces using split,join trick
    region = ''.join(string.split())
    if region == '':
        return None
    # remove undesired characters that may be present and
    # replace - by :
    # N.B., the syntax for translate() differs between python 2 and 3
    try:
        region = region.translate(None, ",;|!{}()").replace("-", ":")
    except:
        region = region.translate({ord(i): None for i in ",;|!{}()"})
    if len(region) == 0:
        print("oh no!")
        raise argparse.ArgumentTypeError(
            "{} is not a valid region".format(string))
    return region


def writableFile(string):
    """
    Simple function that tests if a given path is writable
    """
    try:
        open(string, 'w').close()
        os.remove(string)
    except:
        msg = "{} file can't be opened for writing".format(string)
        raise argparse.ArgumentTypeError(msg)
    return string


"""
Arguments used by heatmapper and profiler
"""


def heatmapperMatrixArgs(args=None):
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--matrixFile', '-m',
                          help='Matrix file from the computeMatrix tool.',
                          type=argparse.FileType('r'),
                          )

    required.add_argument('--outFileName', '-out',
                          help='File name to save the image to. The file '
                          'ending will be used to determine the image '
                          'format. The available options are: "png", '
                          '"eps", "pdf" and "svg", e.g., MyHeatmap.png.',
                          type=writableFile,
                          required=True)
    return parser


def heatmapperOutputArgs(args=None,
                         mode=['heatmap', 'profile'][0]):
    parser = argparse.ArgumentParser(add_help=False)
    output = parser.add_argument_group('Output options')

    output.add_argument(
        '--outFileSortedRegions',
        help='File name into which the regions are saved '
        'after skipping zeros or min/max threshold values. The '
        'order of the regions in the file follows the sorting '
        'order selected. This is useful, for example, to '
        'generate other heatmaps while keeping the sorting of the '
        'first heatmap. Example: Heatmap1sortedRegions.bed',
        metavar='FILE',
        type=argparse.FileType('w'))

    if mode == 'heatmap':
        output.add_argument('--outFileNameMatrix',
                            help='If this option is given, then the matrix '
                            'of values underlying the heatmap will be saved '
                            'using this name, e.g. MyMatrix.tab.',
                            metavar='FILE',
                            type=writableFile)
    elif mode == 'profile':
        output.add_argument('--outFileNameData',
                            help='File name to save the data '
                            'underlying data for the average profile, e.g. '
                            'myProfile.tab.',
                            type=writableFile)
    output.add_argument(
        '--dpi',
        help='Set the DPI to save the figure.',
        type=int,
        default=200)

    return parser


def heatmapperOptionalArgs(mode=['heatmap', 'profile'][0]):

    parser = argparse.ArgumentParser(add_help=False)
    cluster = parser.add_argument_group('Clustering arguments')
    cluster.add_argument(
        '--kmeans',
        help='Number of clusters to compute. When this '
        'option is set, the matrix is split into clusters '
        'using the k-means algorithm. Only works for data that '
        'is not grouped, otherwise only the first group will '
        'be clustered. If more specific clustering methods '
        'are required, then save the underlying matrix '
        'and run the clustering using other software. The plotting  '
        'of the clustering may fail with an error if a '
        'cluster has very few members compared to the total number '
        'or regions.',
        type=int)
    cluster.add_argument(
        '--hclust',
        help='Number of clusters to compute. When this '
        'option is set, then the matrix is split into clusters '
        'using the hierarchical clustering algorithm, using "ward linkage". '
        'Only works for data that is not grouped, otherwise only the first '
        'group will be clustered. --hclust could be very slow if you have '
        '>1000 regions. In those cases, you might prefer --kmeans or if more '
        'clustering methods are required you can save the underlying matrix and run '
        'the clustering using  other software. The plotting of the clustering may '
        'fail with an error if a cluster has very few members compared to the '
        'total number of regions.',
        type=int)

    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument("--help", "-h", action="help",
                          help="show this help message and exit")
    optional.add_argument('--version', action='version',
                          version='%(prog)s {}'.format(__version__))
    if mode == 'profile':
        optional.add_argument(
            '--averageType',
            default='mean',
            choices=["mean", "median", "min",
                     "max", "std", "sum"],
            help='The type of statistic that should be used for the '
            'profile. The options are: "mean", "median", "min", "max", '
            '"sum" and "std".')

        optional.add_argument('--plotHeight',
                              help='Plot height in cm.',
                              type=float,
                              default=7)

        optional.add_argument('--plotWidth',
                              help='Plot width in cm. The minimum value is 1 cm.',
                              type=float,
                              default=11)

        optional.add_argument(
            '--plotType',
            help='"lines" will plot the profile line based '
            'on the average type selected. "fill" '
            'fills the region between zero and the profile '
            'curve. The fill in color is semi transparent to '
            'distinguish different profiles. "se" and "std" '
            'color the region between the profile and the '
            'standard error or standard deviation of the data. '
            'As in the case of '
            'fill, a semi-transparent color is used. '
            '"overlapped_lines" plots each region\'s value, one on '
            'top of the other. "heatmap" plots a '
            'summary heatmap.',
            choices=['lines', 'fill', 'se', 'std', 'overlapped_lines', 'heatmap'],
            default='lines')

        optional.add_argument('--colors',
                              help='List of colors to use '
                              'for the plotted lines (N.B., not applicable to \'--plotType overlapped_lines\'). Color names '
                              'and html hex strings (e.g., #eeff22) '
                              'are accepted. The color names should '
                              'be space separated. For example, '
                              '--colors red blue green ',
                              nargs='+')

        optional.add_argument('--numPlotsPerRow',
                              help='Number of plots per row',
                              type=int,
                              default=8)

    elif mode == 'heatmap':
        optional.add_argument('--sortRegions',
                              help='Whether the heatmap should present '
                              'the regions sorted. The default is '
                              'to sort in descending order based on '
                              'the mean value per region.',
                              choices=["descend", "ascend", "no"],
                              default='descend')

        optional.add_argument('--sortUsing',
                              help='Indicate which method should be used for '
                              'sorting. For each row the '
                              'method is computed.',
                              choices=["mean", "median", "max", "min", "sum",
                                       "region_length"],
                              default='mean')

        optional.add_argument('--sortUsingSamples',
                              help='List of sample numbers (order as in matrix), '
                              'that are used for sorting by --sortUsing, '
                              'no value uses all samples, '
                              'example: --sortUsingSamples 1 3',
                              type=int, nargs='+')

        optional.add_argument(
            '--averageTypeSummaryPlot',
            default='mean',
            choices=["mean", "median", "min",
                     "max", "std", "sum"],
            help='Define the type of statistic that should be plotted in the '
            'summary image above the heatmap. The options are: "mean", '
            '"median", "min", "max", "sum" and "std".')

        optional.add_argument(
            '--missingDataColor',
            default='black',
            help='If --missingDataAsZero is not set, such cases '
            'will be colored in black by default. Using this '
            'parameter, a different color can be set. A value '
            'between 0 and 1 will be used for a gray scale '
            '(black is 0). For a list of possible color '
            'names see: http://packages.python.org/ete2/'
            'reference/reference_svgcolors.html. '
            'Other colors can be specified using the #rrggbb '
            'notation.')

        from matplotlib import cm
        color_options = "', '".join([m for m in cm.datad
                                     if not m.endswith('_r')])

        optional.add_argument(
            '--colorMap',
            help='Color map to use for the heatmap. If more than one heatmap is being plotted the color '
                 'of each heatmap can be enter individually (e.g. `--colorMap Reds Blues`). Color maps '
                 'are recycled if the number of color maps is smaller than the number of heatmaps being '
                 'plotted. Available values can be seen here: http://matplotlib.org/users/colormaps.html '
                 'The available options are: \'' + color_options + '\'',
            default=['RdYlBu'],
            nargs='+')

        optional.add_argument(
            '--alpha',
            default=1.0,
            type=check_float_0_1,
            help='The alpha channel (transparency) to use for the heatmaps. The default is 1.0 and values '
                 'must be between 0 and 1.')

        optional.add_argument(
            '--colorList',
            help='List of colors to use to create a colormap. For example, if `--colorList black,yellow,blue` '
                 'is set (colors separated by comas) then a color map that starts with black, continues to '
                 'yellow and finishes in blue is created. If this option is selected, it overrides the --colorMap '
                 'chosen. The list of valid color names can be seen here: '
                 'http://matplotlib.org/examples/color/named_colors.html  '
                 'Hex colors are valid (e.g #34a2b1). If individual colors for different heatmaps '
                 'need to be specified they need to be separated by space as for example: '
                 '`--colorList "white,#cccccc" "white,darkred"` '
                 'As for --colorMap, the color lists are recycled if their number is smaller thatn the number of'
                 'plotted heatmaps.  '
                 'The number of transitions is defined by the --colorNumber option.',

            nargs='+')

        optional.add_argument(
            '--colorNumber',
            help='N.B., --colorList is required for an effect. This controls the '
            'number of transitions from one color to the other. If --colorNumber is '
            'the number of colors in --colorList then there will be no transitions '
            'between the colors.',
            type=int,
            default=256)

        optional.add_argument('--zMin', '-min',
                              default=None,
                              help='Minimum value for the heatmap intensities. Multiple values, separated by '
                                   'spaces can be set for each heatmap. If the number of zMin values is smaller than'
                                   'the number of heatmaps the values are recycled.',
                              type=float,
                              nargs='+')
        optional.add_argument('--zMax', '-max',
                              default=None,
                              help='Maximum value for the heatmap intensities. Multiple values, separated by '
                                   'spaces can be set for each heatmap. If the number of zMax values is smaller than'
                                   'the number of heatmaps the values are recycled.',
                              type=float,
                              nargs='+')
        optional.add_argument('--heatmapHeight',
                              help='Plot height in cm. The default for the heatmap '
                              'height is 28. The minimum value is '
                              '3 and the maximum is 100.',
                              type=float,
                              default=28)

        optional.add_argument('--heatmapWidth',
                              help='Plot width in cm. The default value is 4 '
                              'The minimum value is 1 and the '
                              'maximum is 100.',
                              type=float,
                              default=4)
        optional.add_argument(
            '--whatToShow',
            help='The default is to include a summary or profile plot on top '
            'of the heatmap and a heatmap colorbar. Other options are: '
            '"plot and heatmap", "heatmap only", "heatmap and '
            'colorbar", and the default "plot, heatmap and '
            'colorbar".',
            choices=["plot, heatmap and colorbar",
                     "plot and heatmap", "heatmap only",
                     "heatmap and colorbar"],
            default='plot, heatmap and colorbar')

        optional.add_argument(
            '--boxAroundHeatmaps',
            help='By default black boxes are plot around heatmaps. This can be turned off '
                 'by setting --boxAroundHeatmaps no',
            default='yes')

        optional.add_argument('--xAxisLabel', '-x',
                              default='gene distance (bp)',
                              help='Description for the x-axis label.')

    # end elif
    optional.add_argument('--startLabel',
                          default='TSS',
                          help='[only for scale-regions mode] Label shown '
                          'in the plot for the start of '
                          'the region. Default is TSS (transcription '
                          'start site), but could be changed to anything, '
                          'e.g. "peak start". '
                          'Same for the --endLabel option. See below.')
    optional.add_argument('--endLabel',
                          default='TES',
                          help='[only for scale-regions mode] Label '
                          'shown in the plot for the region '
                          'end. Default is TES (transcription end site).')
    optional.add_argument('--refPointLabel',
                          help='[only for reference-point mode] Label '
                          'shown in the plot for the '
                          'reference-point. Default '
                          'is the same as the reference point selected '
                          '(e.g. TSS), but could be anything, e.g. '
                          '"peak start".',
                          default='TSS')
    optional.add_argument('--nanAfterEnd',
                          help=argparse.SUPPRESS,
                          default=False)

    optional.add_argument('--regionsLabel', '-z',
                          help='Labels for the regions plotted in the '
                          'heatmap. If more than one region is being '
                          'plotted, a list of labels separated by spaces is required. '
                          'If a label itself contains a space, then quotes are '
                          'needed. For example, --regionsLabel label_1, "label 2". ',
                          nargs='+')

    optional.add_argument('--samplesLabel',
                          help='Labels for the samples plotted. The '
                          'default is to use the file name of the '
                          'sample. The sample labels should be separated '
                          'by spaces and quoted if a label itself'
                          'contains a space E.g. --samplesLabel label-1 "label 2"  ',
                          nargs='+')

    optional.add_argument('--plotTitle', '-T',
                          help='Title of the plot, to be printed on top of '
                          'the generated image. Leave blank for no title.',
                          default='')

    optional.add_argument('--yAxisLabel', '-y',
                          default='',
                          help='Y-axis label for the top panel.')

    optional.add_argument('--yMin',
                          default=None,
                          help='Minimum value for the Y-axis.')
    optional.add_argument('--yMax',
                          default=None,
                          help='Maximum value for the Y-axis.')

    optional.add_argument('--legendLocation',
                          default='best',
                          choices=['best',
                                   'upper-right',
                                   'upper-left',
                                   'upper-center',
                                   'lower-left',
                                   'lower-right',
                                   'lower-center',
                                   'center',
                                   'center-left',
                                   'center-right',
                                   'none'
                                   ],
                          help='Location for the legend in the summary plot. '
                               'Note that "none" does not work for the profiler.')

    optional.add_argument('--perGroup',
                          help='The default is to plot all groups of regions by '
                          'sample. Using this option instead plots all samples by '
                          'group of regions. Note that this is only useful if you '
                          'have multiple groups of regions. by sample rather than '
                          'group.',
                          action='store_true')

    optional.add_argument('--plotFileFormat',
                          metavar='',
                          help='Image format type. If given, this '
                          'option overrides the '
                          'image format based on the plotFile ending. '
                          'The available options are: "png", '
                          '"eps", "pdf" and "svg"',
                          choices=['png', 'pdf', 'svg', 'eps'])

    optional.add_argument('--verbose',
                          help='If set, warning messages and '
                          'additional information are given.',
                          action='store_true')
    return parser
