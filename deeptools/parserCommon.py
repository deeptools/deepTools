import sys
import argparse
import deeptools.config as cfg
from deeptools._version import __version__


def output(args=None):
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Output')
    group.add_argument('--outFileName', '-o',
                       help='Output file name.',
                       metavar='FILENAME',
                       type=writableFile,
                       required=True)

    default = checkBigWig('bigwig')
    group.add_argument('--outFileFormat', '-of',
                       help='Output file type. Either "bigwig" or "bedgraph".',
                       choices=['bigwig', 'bedgraph'],
                       default=default)

    return parser


def read_options():
    """Common arguments related to bam files and the interpretation
    of the read coverage
    """
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group('Read processing options')

    group.add_argument('--fragmentLength', '-f',
                       help='Only when --extendPairedEnds is set. '
                            'Length of the average fragment size. Reads will '
                            'be extended to match this length unless they are '
                            'paired-end, in which case they will be extended to '
                            'match the fragment length. If this value is set to '
                            'between 1 and the read length, the read will not be '
                            'extended. *NOTE*: If the BAM files contain mated and '
                            'unmated paired-end reads, unmated reads will be '
                            'extended to match the --fragmentLength.',
                       type=int,
                       metavar="INT bp")

    group.add_argument('--extendPairedEnds',
                       help='If set, paired end reads are extended to match the '
                            'fragment length between the mate reads. By default '
                            '*each* read mate is extended. This can be modified using'
                            'the sam flags (see samFlagInclude and samFlagExclude '
                            'options) to keep only the first or the second mate. Mate '
                            'reads that are in different chromosomes or that are to far'
                            'apart are extended to the given --fragmentLength',
                       action='store_true')

    group.add_argument('--ignoreDuplicates',
                       help='If set, reads that have the same orientation '
                       'and start position will be considered only '
                       'once. If reads are paired, the mate position '
                       'also has to coincide to ignore a read.',
                       action='store_true'
                       )

    group.add_argument('--minMappingQuality',
                       metavar='INT',
                       help='If set, only reads that have a mapping '
                       'quality score equal or higher than --minMappingQuality are '
                       'considered.',
                       type=int,
                       )

    group.add_argument('--centerReads',
                       help='By adding this option reads are centered with '
                       'respect to the fragment length. For paired-end data '
                       'the read is centered at the fragment length defined '
                       'by the two fragment ends. For single-end data, the '
                       'given fragment length is used. This option is '
                       'useful to get a sharper signal around enriched '
                       'regions.',
                       action='store_true')

    group.add_argument('--samFlagInclude',
                       help='Include reads based on the SAM flag. For example, '
                       'to get only reads that are the first mate use a flag of 64. '
                       'This is useful to count properly paired reads only once, '
                       'otherwise the second mate will be also considered for the '
                       'coverage.',
                       metavar= 'INT',
                       default=None,
                       type=int,
                       required=False)

    group.add_argument('--samFlagExclude',
                       help='Exclude reads based on the SAM flag. For example, '
                       'to get only reads that map to the forward strand, use '
                       '--samFlagExclude 16. Where 16 is the sam flag for reads '
                       'that map to the reverse strand.',
                       metavar= 'INT',
                       default=None,
                       type=int,
                       required=False)

    return parser


def getParentArgParse(args=None, binSize=True):
    """
    Typical arguments for several tools
    """

    parser = argparse.ArgumentParser(add_help=False)
    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument('--version', action='version',
                          version='%(prog)s {}'.format(__version__))

    if binSize:
        optional.add_argument('--binSize', '-bs',
                              help='Size of the bins in bp for the output '
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
                                             "computer is {}.".format(string,e,
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
    region = region.translate(None, ",;|!{}()").replace("-", ":")
    if len(region) == 0:
        raise argparse.ArgumentTypeError(
            "{} is not a valid region".format(string))
    return region


def writableFile(string):
    """
    Simple function that tests if a given path is writable
    """
    try:
        open(string, 'w').close()
    except:
        msg = "{} file can be opened for writing".format(string)
        raise argparse.ArgumentTypeError(msg)
    return string


def checkBigWig(string):
    """
    Checks if the path to USCS bedGraphToBigWig as set in the config
    is installed and is executable.
    """
    if string == 'bigwig':
        bedgraph_to_bigwig = cfg.config.get('external_tools',
                                            'bedgraph_to_bigwig')
        if not cfg.checkProgram(bedgraph_to_bigwig, 'h',
                                'http://hgdownload.cse.ucsc.edu/admin/exe/'):
            msg = "The output is set by default to 'bedgraph'\n"
            sys.stderr.write(msg)
            return 'bedgraph'

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
                        help='File name to save the image. The file '
                        'ending will be used to determine the image '
                        'format. The available options are: "png", "emf", '
                        '"eps", "pdf" and "svg", e. g. MyHeatmap.png.',
                        type=writableFile,
                        required=True)
    return parser


def heatmapperOutputArgs(args=None,
                         mode=['heatmap', 'profile'][0]):
    parser = argparse.ArgumentParser(add_help=False)
    output = parser.add_argument_group('Output options')

    output.add_argument('--outFileNameData',
                        help='File name to save the data '
                        'underlying data for the average profile, e.g. '
                        'myProfile.tab.',
                        type=argparse.FileType('w'))

    output.add_argument(
        '--outFileSortedRegions',
        help='File name in which the regions are saved '
        'after skipping zeros or min/max threshold values. The '
        'order of the regions in the file follows the sorting '
        'order selected. This is useful, for example, to '
        'generate other heatmaps keeping the sorting of the '
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
    return parser


def heatmapperOptionalArgs(mode=['heatmap', 'profile'][0]):

    parser = argparse.ArgumentParser(add_help=False)
    cluster = parser.add_argument_group('Clustering arguments')
    cluster.add_argument(
        '--kmeans',
        help='Number of clusters to compute. When this '
        'option is set, then the matrix is split into clusters '
        'using the kmeans algorithm. Only works for data that '
        'is not grouped, otherwise only the first group will '
        'be clustered. If more specific clustering methods '
        'are required it is advisable to save the underlying matrix '
        'and run the clustering using other software. The plotting  '
        'of the clustering may fail (Error: Segmentation fault) if a '
        'cluster has very few members compared to the total number '
        'or regions.',
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
            help='Define the type of statistic that should be used for the '
            'profile. The options are: "mean", "median", "min", "max", '
            '"sum" and "std".')

        optional.add_argument('--plotHeight',
                              help='height in cm.',
                              type=float,
                              default=7)

        optional.add_argument('--plotWidth',
                              help='Width in cm.  The minimum value is 1 cm.',
                              type=float,
                              default=11)

        optional.add_argument(
            '--plotType',
            help='For the summary plot (profile) only. The '
            '"lines" option will plot the profile line based '
            'on the average type selected. The "fill" option '
            'fills the region between zero and the profile '
            'curve. The fill in color is semi transparent to '
            'distinguish different profiles. The "se" and "std" options '
            'colors the region between the profile and the '
            'standard error or standard deviation of the data. '
            'As in the case of '
            'fill, a semi-transparent color is used. The option '
            '"overlapped_lines" plots each region values, one on '
            'top of the other; this option only works if '
            '--onePlotPerGroup is set. The option "heatmap" plots a '
            'summary heatmap.',
            choices=['lines', 'fill', 'se', 'std', 'overlapped_lines', 'heatmap'],
            default='lines')

        optional.add_argument('--colors',
                              help='List of colors to use '
                              'for the plotted lines. Color names '
                              'and html hex strings (e.g. #eeff22) '
                              'are accepted. The color names should '
                              'be given separated by spaces. For example '
                              '--colors red blue green ',
                              nargs='+')

        optional.add_argument('--numPlotsPerRow',
                              help='Number of plots to plot in a row',
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
            'parameter a different color can be set. A value '
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
            '--colorMap', default='RdYlBu',
            help='Color map to use for the heatmap. Available values can be '
            'seen here: '
            'http://www.astro.lsa.umich.edu/~msshin/science/code/'
            'matplotlib_cm/ The available options are: \'' +
            color_options + '\'')

        optional.add_argument(
            '--colorList',
            help='List of color to create a colormap. For example if  '
            '--colorList black yellow blue is set (colors separated by '
            'spaces) then a color map that starts with black, continues '
            'to yellow and finishes in blue is created. If this option is'
            'selected, it overrides the --colorMap selected.',
            nargs='+')

        optional.add_argument(
            '--colorNumber',
            help='if --colorList is given, colorNumber controls the number of '
            'transitions from one color to the other. If --colorNumber is set '
            'to the same number of colors as in --colorList, then per each '
            'color a discrete interval  is shown (i.e. no transitions)',
            type=int,
            default=256)

        optional.add_argument('--zMin', '-min',
                              default=None,
                              help='Minimum value for the heatmap intensities.')
        optional.add_argument('--zMax', '-max',
                              default=None,
                              help='Maximum value for the heatmap intensities.')
        optional.add_argument('--heatmapHeight',
                              help='Height in cm. The default for the heatmap '
                              'height is 25. The minimum value is '
                              '3 and the maximum is 100.',
                              type=float,
                              default=28)

        optional.add_argument('--heatmapWidth',
                              help='Width in cm. The default value is 4 '
                              'The minimum value is 1 and the '
                              'maximum is 100.',
                              type=float,
                              default=4)
        optional.add_argument(
            '--whatToShow',
            help='The default is to include a summary or profile plot on top '
            'of the heatmap and a heatmap colorbar. Other options are: '
            '"plot and heatmap", '
            '"heatmap only", "colorbar only", "heatmap and '
            'colorbar", and the default "plot, heatmap and '
            'colorbar".',
            choices=["plot, heatmap and colorbar",
                     "plot and heatmap", "heatmap only",
                     "colorbar only", "heatmap and colorbar"],
            default='plot, heatmap and colorbar')

    ## end elif
    optional.add_argument('--startLabel',
                          default='TSS',
                          help='[only for scale-regions mode] Label shown '
                          'in the plot for the start of '
                          'the region. Default is TSS (transcription '
                          'start site), but could be changed to anything, '
                          'e.g. "peak start".'
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
                          default='genes',
                          help='Labels for the regions plotted in the '
                          'heatmap. If more than one region is being '
                          'plotted a list of lables '
                          'separated by comma and limited by quotes, '
                          'is requires. For example, '
                          ' --regionsLabel "label1, label2". '
                          'Default is "genes".')

    optional.add_argument('--samplesLabel',
                          help='Labels for the samples plotted. The '
                          'default is to use the file name of the '
                          'sample. The sample names should be separated '
                               'by spaces and quoted if they label itself'
                               'contains a space E.g. --samplesLabel label-1 "label 2"  ',
                          nargs='+')

    optional.add_argument('--plotTitle', '-T',
                          help='Title of the plot, to be printed on top of '
                          'the generated image. Leave blank for no title.',
                          default='')

    optional.add_argument('--xAxisLabel', '-x',
                          default='gene distance (bp)',
                          help='Description for the x-axis label.')
    optional.add_argument('--yAxisLabel', '-y',
                          default='',
                          help='Description for the y-axis label for the '
                          'top panel.')

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
                          help='Location for the legend in the summary plot')

    optional.add_argument('--perGroup',
                          help='The default is to combine all group '
                          'plots into one. If multiple samples are '
                          'present, then for each sample a new plot is '
                          'made. If --perGroup is set, then, for each '
                          'group, all data from different samples '
                          'is combined in to one plot.',
                          action='store_true')

    optional.add_argument('--plotFileFormat',
                          metavar='',
                          help='image format type. If given, this '
                          'option overrides the '
                          'image format based on the plotFile ending. '
                          'The available options are: "png", "emf", '
                          '"eps", "pdf" and "svg"',
                          choices=['png', 'pdf', 'svg', 'eps', 'emf'])

    optional.add_argument('--verbose',
                          help='If set warning messages and '
                          'addition information are given.',
                          action='store_true')
    return parser


def bam_total_reads(bam_handle, chroms_to_ignore):
    """Count the total number of mapped reads in a BAM file, filtering
    the chromosome given in chroms_to_ignore list
    """
    if chroms_to_ignore:
        import pysam

        lines = pysam.idxstats(bam_handle.filename)
        tot_mapped_reads = 0
        for line in lines:
            chrom, _len, nmapped, _nunmapped = line.split('\t')
            if chrom not in chroms_to_ignore:
                tot_mapped_reads += int(nmapped)

    else:
        tot_mapped_reads = bam_handle.mapped

    return tot_mapped_reads
