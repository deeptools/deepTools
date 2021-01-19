#!/usr/bin/env python
import deeptools.heatmapper as heatmapper
import deeptoolsintervals.parse as dti
import numpy as np
import argparse
import sys
import os
import csv
from deeptools._version import __version__


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
This tool performs a variety of operations on files produced by computeMatrix.

detailed help:

  computeMatrixOperations info -h

or

  computeMatrixOperations relabel -h

or

  computeMatrixOperations subset -h

or

  computeMatrixOperations filterStrand -h

or

  computeMatrixOperations filterValues -h

or

  computeMatrixOperations rbind -h

or

  computeMatrixOperations cbind -h

or
  computeMatrixOperations sort -h

or
  computeMatrixOperations dataRange -h

""",
        epilog='example usages:\n'
               'computeMatrixOperations subset -m input.mat.gz -o output.mat.gz --group "group 1" "group 2" --samples "sample 3" "sample 10"\n\n'
               ' \n\n')

    subparsers = parser.add_subparsers(
        title='Commands',
        dest='command',
        metavar='')

    # info
    subparsers.add_parser(
        'info',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[infoArgs()],
        help="Print group and sample information",
        usage='An example usage is:\n  computeMatrixOperations info -m input.mat.gz\n\n')

    # relabel
    subparsers.add_parser(
        'relabel',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[infoArgs(), relabelArgs()],
        help="Change sample and/or group label information",
        usage='An example usage is:\n  computeMatrixOperations relabel -m input.mat.gz -o output.mat.gz --sampleLabels "sample 1" "sample 2"\n\n')

    # subset
    subparsers.add_parser(
        'subset',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[infoArgs(), subsetArgs()],
        help="Actually subset the matrix. The group and sample orders are honored, so one can also reorder files.",
        usage='An example usage is:\n  computeMatrixOperations subset -m '
        'input.mat.gz -o output.mat.gz --groups "group 1" "group 2" '
        '--samples "sample 3" "sample 10"\n\n')

    # filterStrand
    subparsers.add_parser(
        'filterStrand',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[infoArgs(), filterStrandArgs()],
        help="Filter entries by strand.",
        usage='Example usage:\n  computeMatrixOperations filterStrand -m '
        'input.mat.gz -o output.mat.gz --strand +\n\n')

    # filterValues
    subparsers.add_parser(
        'filterValues',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[infoArgs(), filterValuesArgs()],
        help="Filter entries by min/max value.",
        usage='Example usage:\n  computeMatrixOperations filterValues -m '
        'input.mat.gz -o output.mat.gz --min 10 --max 1000\n\n')

    # rbind
    subparsers.add_parser(
        'rbind',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[bindArgs()],
        help="merge multiple matrices by concatenating them head to tail. This assumes that the same samples are present in each in the same order.",
        usage='Example usage:\n  computeMatrixOperations rbind -m '
        'input1.mat.gz input2.mat.gz -o output.mat.gz\n\n')

    # cbind
    subparsers.add_parser(
        'cbind',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[bindArgs()],
        help="merge multiple matrices by concatenating them left to right. No assumptions are made about the row order. Regions not present in the first file specified are ignored. Regions missing in subsequent files will result in NAs. Regions are matches based on the first 6 columns of the computeMatrix output (essentially the columns in a BED file).",
        usage='Example usage:\n  computeMatrixOperations cbind -m '
        'input1.mat.gz input2.mat.gz -o output.mat.gz\n\n')

    # sort
    subparsers.add_parser(
        'sort',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[sortArgs()],
        help='Sort a matrix file to correspond to the order of entries in the desired input file(s). The groups of regions designated by the files must be present in the order found in the output of computeMatrix (otherwise, use the subset command first). Note that this subcommand can also be used to remove unwanted regions, since regions not present in the input file(s) will be omitted from the output.',
        usage='Example usage:\n  computeMatrixOperations sort -m input.mat.gz -R regions1.bed regions2.bed regions3.gtf -o input.sorted.mat.gz\n\n')

    # dataRange
    subparsers.add_parser(
        'dataRange',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[infoArgs()],
        help='Returns the min, max, median, 10th and 90th percentile of the matrix values per sample.',
        usage='Example usage:\n  computeMatrixOperations dataRange -m input.mat.gz\n\n')

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def bindArgs():
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    required.add_argument('--matrixFile', '-m',
                          help='Matrix files from the computeMatrix tool.',
                          nargs='+',
                          required=True)

    required.add_argument('--outFileName', '-o',
                          help='Output file name',
                          required=True)

    return parser


def infoArgs():
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    required.add_argument('--matrixFile', '-m',
                          help='Matrix file from the computeMatrix tool.',
                          required=True)

    return parser


def relabelArgs():
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    required.add_argument('--outFileName', '-o',
                          help='Output file name',
                          required=True)

    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument('--groupLabels',
                          nargs='+',
                          help="Groups labels. If none are specified then the current labels will be kept.")

    optional.add_argument('--sampleLabels',
                          nargs='+',
                          help="Sample labels. If none are specified then the current labels will be kept.")

    return parser


def subsetArgs():
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    required.add_argument('--outFileName', '-o',
                          help='Output file name',
                          required=True)

    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument('--groups',
                          nargs='+',
                          help="Groups to include. If none are specified then all will be included.")

    optional.add_argument('--samples',
                          nargs='+',
                          help="Samples to include. If none are specified then all will be included.")

    return parser


def filterStrandArgs():
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    required.add_argument('--outFileName', '-o',
                          help='Output file name',
                          required=True)

    required.add_argument('--strand', '-s',
                          help='Strand',
                          choices=['+', '-', '.'],
                          required=True)

    return parser


def filterValuesArgs():
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    required.add_argument('--outFileName', '-o',
                          help='Output file name',
                          required=True)

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--min',
                          help='Minimum value. Any row having a single entry less than this will be excluded. The default is no minimum.',
                          type=float,
                          default=None)

    optional.add_argument('--max',
                          help='Maximum value. Any row having a single entry more than this will be excluded. The default is no maximum.',
                          type=float,
                          default=None)

    return parser


def sortArgs():
    parser = argparse.ArgumentParser(add_help=False)
    required = parser.add_argument_group('Required arguments')

    required.add_argument('--matrixFile', '-m',
                          help='Matrix file from the computeMatrix tool.',
                          required=True)

    required.add_argument('--outFileName', '-o',
                          help='Output file name',
                          required=True)

    required.add_argument('--regionsFileName', '-R',
                          help='File name(s), in BED or GTF format, containing the regions. '
                               'If multiple bed files are given, each one is '
                               'considered a group that can be plotted separately. '
                               'Also, adding a "#" symbol in the bed file causes all '
                               'the regions until the previous "#" to be considered '
                               'one group. Alternatively for BED files, putting '
                               'deepTools_group in the header can be used to indicate a '
                               'column with group labels. Note that these should be '
                               'sorted such that all group entries are together.',
                          required=True,
                          nargs='+')

    optional = parser.add_argument_group('Optional arguments')

    optional.add_argument('--transcriptID',
                          default='transcript',
                          help='When a GTF file is used to provide regions, only '
                          'entries with this value as their feature (column 3) '
                          'will be processed as transcripts. (Default: %(default)s)')

    optional.add_argument('--transcript_id_designator',
                          default='transcript_id',
                          help='Each region has an ID (e.g., ACTB) assigned to it, '
                          'which for BED files is either column 4 (if it exists) '
                          'or the interval bounds. For GTF files this is instead '
                          'stored in the last column as a key:value pair (e.g., as '
                          '\'transcript_id "ACTB"\', for a key of transcript_id '
                          'and a value of ACTB). In some cases it can be '
                          'convenient to use a different identifier. To do so, set '
                          'this to the desired key. (Default: %(default)s)')

    return parser


def printInfo(matrix):
    """
    Print the groups and samples
    """

    print("Groups:")
    for group in matrix.matrix.group_labels:
        print("\t{0}".format(group))

    print("Samples:")
    for sample in matrix.matrix.sample_labels:
        print("\t{0}".format(sample))


def printDataRange(matrix):
    """
    Prints the min, max, median, 10th and 90th percentile of the matrix values per sample.
    """
    print("Samples\tMin\tMax\tMedian\t10th\t90th")
    for i, sample in enumerate(matrix.matrix.sample_labels):
        start = matrix.matrix.sample_boundaries[i]
        end = matrix.matrix.sample_boundaries[i + 1]
        sample_matrix = matrix.matrix.matrix[..., start:end]
        print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(sample, np.amin(sample_matrix),
                                                    np.amax(sample_matrix),
                                                    np.ma.median(sample_matrix),
                                                    np.percentile(sample_matrix, 10),
                                                    np.percentile(sample_matrix, 90)))


def relabelMatrix(matrix, args):
    """
    Relabel the samples and groups in a matrix
    """
    if args.groupLabels:
        if len(args.groupLabels) != len(matrix.matrix.group_labels):
            sys.exit("You specified {} group labels, but {} are required.\n".format(len(args.groupLabels), len(matrix.matrix.group_labels)))
        matrix.matrix.group_labels = args.groupLabels
    if args.sampleLabels:
        if len(args.sampleLabels) != len(matrix.matrix.sample_labels):
            sys.exit("You specified {} sample labels, but {} are required.\n".format(len(args.sampleLabels), len(matrix.matrix.sample_labels)))
        matrix.matrix.sample_labels = args.sampleLabels


def getGroupBounds(args, matrix):
    """
    Given the group labels, return an indexing array and the resulting boundaries
    """
    bounds = matrix.parameters['group_boundaries']
    if args.groups is None:
        return range(0, matrix.matrix.matrix.shape[0]), np.array(bounds)
    else:
        o = list()
        obounds = [0]
        for group in args.groups:
            if group not in matrix.matrix.group_labels:
                sys.exit("Error: '{0}' is not a valid group\n".format(group))
            idx = matrix.matrix.group_labels.index(group)
            o.extend(range(bounds[idx], bounds[idx + 1]))
            obounds.append(bounds[idx + 1] - bounds[idx])
        return o, np.cumsum(obounds)


def getSampleBounds(args, matrix):
    """
    Given the sample labels, return an indexing array
    """
    bounds = matrix.parameters['sample_boundaries']
    if args.samples is None:
        return np.arange(0, matrix.matrix.matrix.shape[1])
    else:
        o = list()
        for sample in args.samples:
            if sample not in matrix.matrix.sample_labels:
                sys.exit("Error: '{0}' is not a valid sample\n".format(sample))
            idx = matrix.matrix.sample_labels.index(sample)
            o.extend(range(bounds[idx], bounds[idx + 1]))
        return o


def subsetRegions(hm, bounds):
    out = []
    for x in bounds:
        reg = hm.matrix.regions[x]
        # we need to add a list of [chrom, [(start, end), (start, end)], name, 0, strand, score)]
        if isinstance(reg, dict):
            # This happens on occasion
            starts = reg["start"].split(",")
            starts = [int(x) for x in starts]
            ends = reg["end"].split(",")
            ends = [int(x) for x in ends]
            regs = [(x, y) for x, y in zip(starts, ends)]
            out.append([reg["chrom"], regs, reg["name"], 0, reg["strand"], reg["score"]])
        else:
            out.append(reg)
    return out


def filterHeatmap(hm, args):
    bounds = [0]
    regions = []
    keep = []
    for region in hm.matrix.regions:
        if region[4] == args.strand:
            keep.append(True)
            regions.append(region)
        else:
            keep.append(False)
    keep = np.array(keep)

    # Get the new bounds
    for idx in range(1, len(hm.matrix.group_boundaries)):
        i = int(np.sum(keep[hm.matrix.group_boundaries[idx - 1]:hm.matrix.group_boundaries[idx]]))
        bounds.append(bounds[idx - 1] + i)

    hm.matrix.group_boundaries = bounds

    # subset the matrix
    hm.matrix.matrix = hm.matrix.matrix[keep, :]
    hm.matrix.regions = regions


def filterHeatmapValues(hm, minVal, maxVal):
    bounds = [0]
    regions = []
    keep = []
    if minVal is None:
        minVal = -np.inf
    if maxVal is None:
        maxVal = np.inf
    np.warnings.filterwarnings('ignore')
    for i, (x, y) in enumerate(zip(np.nanmin(hm.matrix.matrix, axis=1), np.nanmax(hm.matrix.matrix, axis=1))):
        # x/y will be nan iff a row is entirely nan. Don't filter.
        if np.isnan(x) or (x >= minVal and y <= maxVal):
            keep.append(True)
            regions.append(hm.matrix.regions[i])
        else:
            keep.append(False)
    keep = np.array(keep)

    # Get the new bounds
    for idx in range(1, len(hm.matrix.group_boundaries)):
        i = int(np.sum(keep[hm.matrix.group_boundaries[idx - 1]:hm.matrix.group_boundaries[idx]]))
        bounds.append(bounds[idx - 1] + i)

    hm.matrix.group_boundaries = bounds

    # subset the matrix
    hm.matrix.matrix = hm.matrix.matrix[keep, :]
    hm.matrix.regions = regions


def insertMatrix(hm, hm2, groupName):
    """
    Given two heatmapper objects and a region group name, insert the regions and
    values from hm2 for that group to the end of those for hm.
    """
    # get the bounds for hm
    idx = hm.parameters["group_labels"].index(groupName)
    hmEnd = hm.parameters["group_boundaries"][idx + 1]
    # get the bounds for hm2
    idx2 = hm2.parameters["group_labels"].index(groupName)
    hm2Start = hm2.parameters["group_boundaries"][idx2]
    hm2End = hm2.parameters["group_boundaries"][idx2 + 1]

    # Insert the subset hm2 into hm along axis 0
    hm.matrix.matrix = np.insert(hm.matrix.matrix, hmEnd, hm2.matrix.matrix[hm2Start:hm2End, :], axis=0)

    # Insert the regions
    hm.matrix.regions[hmEnd:hmEnd] = hm2.matrix.regions[hm2Start:hm2End]

    # Increase the group boundaries
    bounds = []
    for idx3, bound in enumerate(hm.parameters["group_boundaries"]):
        if idx3 > idx:
            bound += hm2End - hm2Start
        bounds.append(bound)
    hm.parameters["group_boundaries"] = bounds


def appendMatrix(hm, hm2, groupName):
    """
    Given two heatmapper objects and a region group name, append the values from
    that group in hm2 onto the end of hm.
    """
    # get the bounds for hm2
    idx2 = hm2.parameters["group_labels"].index(groupName)
    hm2Start = hm2.parameters["group_boundaries"][idx2]
    hm2End = hm2.parameters["group_boundaries"][idx2 + 1]

    # Append the matrix
    hm.matrix.matrix = np.concatenate([hm.matrix.matrix, hm2.matrix.matrix[hm2Start:hm2End, :]], axis=0)
    # Update the bounds
    hm.parameters["group_boundaries"].append(hm.parameters["group_boundaries"][-1] + hm2End - hm2Start)
    # Append the regions
    hm.matrix.regions.extend(hm2.matrix.regions[hm2Start:hm2End])


def rbindMatrices(hm, args):
    """
    Bind matrices, top to bottom while accounting for the groups.

    It's assumed that the same samples are present in both and in the exact same order
    """
    hm2 = heatmapper.heatmapper()
    hm.read_matrix_file(args.matrixFile[0])
    for idx in range(1, len(args.matrixFile)):
        hm2.read_matrix_file(args.matrixFile[idx])
        for idx, group in enumerate(hm2.parameters["group_labels"]):
            if group in hm.parameters["group_labels"]:
                insertMatrix(hm, hm2, group)
            else:
                appendMatrix(hm, hm2, group)
                hm.parameters["group_labels"].append(group)

    # Update the group boundaries attribute
    hm.matrix.group_labels = hm.parameters['group_labels']
    hm.matrix.group_boundaries = hm.parameters['group_boundaries']


def cbindMatrices(hm, args):
    """
    Bind columns from different matrices according to the group and region names

    Missing regions are left as NA
    """
    hm2 = heatmapper.heatmapper()

    # Make a dict of region name:row associations
    hm.read_matrix_file(args.matrixFile[0])
    d = dict({x: dict() for x in hm.parameters["group_labels"]})
    for idx, group in enumerate(hm.parameters["group_labels"]):
        s = hm.parameters["group_boundaries"][idx]
        e = hm.parameters["group_boundaries"][idx + 1]
        for idx2, reg in enumerate(hm.matrix.regions[s:e]):
            d[group][reg[2]] = idx2 + s

    # Iterate through the other matrices
    for idx in range(1, len(args.matrixFile)):
        hm2.read_matrix_file(args.matrixFile[idx])
        # Add the sample labels
        hm.parameters['sample_labels'].extend(hm2.parameters['sample_labels'])
        # Add the sample boundaries
        lens = [x + hm.parameters['sample_boundaries'][-1] for x in hm2.parameters['sample_boundaries']][1:]
        hm.parameters['sample_boundaries'].extend(lens)

        # Add on additional NA initialized columns
        ncol = hm.matrix.matrix.shape[1]
        hm.matrix.matrix = np.hstack((hm.matrix.matrix, np.empty(hm2.matrix.matrix.shape)))
        hm.matrix.matrix[:, ncol:] = np.NAN

        # Update the values
        for idx2, group in enumerate(hm2.parameters["group_labels"]):
            if group not in d:
                continue
            s = hm2.parameters["group_boundaries"][idx2]
            e = hm2.parameters["group_boundaries"][idx2 + 1]
            for idx3, reg in enumerate(hm2.matrix.regions[s:e]):
                if reg[2] not in d[group]:
                    continue
                hm.matrix.matrix[d[group][reg[2]], ncol:] = hm2.matrix.matrix[s + idx3, :]

        # Append the special params
        for s in hm.special_params:
            hm.parameters[s].extend(hm2.parameters[s])

    # Update the sample parameters
    hm.matrix.sample_labels = hm.parameters['sample_labels']
    hm.matrix.sample_boundaries = hm.parameters['sample_boundaries']


def loadBED(line, fp, fname, labelColumn, labels, regions, defaultGroup):
    """
    Given a first line, possibly a label column and a list of labels and regions, add the labels and regions in the file to them
    """

    # This is largely parseBED from deeptoolsintervals
    labelIdx = None
    localRegions = {}

    cols = line.strip().split("\t")
    if labelColumn is not None:
        label = cols.pop(labelColumn)
        if label not in labels:
            labels[label] = len(labels)
        labelIdx = labels[label]
        if labelIdx >= len(regions):
            regions.append(localRegions)
        else:
            localRegions = regions[labelIdx]

    if len(cols) >= 6:
        name = cols[3]
    else:
        name = "{0}:{1}-{2}".format(cols[0], cols[1], cols[2])
    localRegions[name] = len(localRegions)

    for line in fp:
        if line.startswith("#") and labelColumn is None:
            if len(localRegions) > 0:
                label = line[1:].strip()
                if len(label):
                    labels[dti.findRandomLabel(labels, label)] = len(labels)
                else:
                    labels[dti.findRandomLabel(labels, os.path.basename(fname))] = len(labels)
                regions.append(localRegions)
                localRegions = dict()
            continue
        elif line.startswith("#") and labelColumn is not None:
            continue

        cols = line.strip().split("\t")
        if len(cols) < 3:
            continue
        if labelColumn is not None:
            label = cols.pop(labelColumn)
            if label not in labels:
                labels[label] = len(labels)
            labelIdx = labels[label]
            if labelIdx >= len(regions):
                regions.append({})
            localRegions = regions[labelIdx]

        if len(cols) >= 6:
            name = cols[3]
        else:
            name = "{0}:{1}-{2}".format(cols[0], cols[1], cols[2])
        name = dti.findRandomLabel(localRegions, name)
        localRegions[name] = len(localRegions)

    # Handle the last group if there is no label
    if labelIdx is None and len(localRegions) > 0:
        if defaultGroup is not None:
            labels[dti.findRandomLabel(labels, defaultGroup)] = len(labels)
        else:
            labels[dti.findRandomLabel(labels, os.path.basename(fname))] = len(labels)
        regions.append(localRegions)


def loadGTFtranscript(cols, label, defaultGroup, transcript_id_designator):
    s = next(csv.reader([cols[8]], delimiter=' '))
    if "deepTools_group" in s and s[-1] != "deepTools_group":
        label = s[s.index("deepTools_group") + 1].rstrip(";")
    elif defaultGroup is not None:
        label = defaultGroup

    if transcript_id_designator not in s or s[-1] == transcript_id_designator:
        sys.stderr.write("Warning: {0} is malformed!\n".format("\t".join(cols)))
        return None, None

    name = s[s.index(transcript_id_designator) + 1].rstrip(";")
    return label, name


def loadGTF(line, fp, fname, labels, regions, transcriptID, transcript_id_designator, defaultGroup):
    """
    Like loadBED, but for a GTF file

    This is largely a copy of what's in deeptoolsintervals
    """
    file_label = dti.findRandomLabel(labels, os.path.basename(fname))

    # handle the first line
    cols = line.split("\t")
    if cols[2].lower() == transcriptID.lower():
        label, name = loadGTFtranscript(cols, file_label, defaultGroup, transcript_id_designator)
        if label is not None:
            if label not in labels:
                labels[label] = len(labels)
                regions.append(dict())
            labelIdx = labels[label]
            regions[labelIdx][name] = len(regions[labelIdx])

    for line in fp:
        if not isinstance(line, str):
            line = line.decode('ascii')
        if not line.startswith('#'):
            cols = line.strip().split('\t')
            if len(cols) == 0:
                continue
            if cols[2].lower() == transcriptID:
                label, name = loadGTFtranscript(cols, file_label, defaultGroup, transcript_id_designator)
                if label is None:
                    continue
                if label not in labels:
                    labels[label] = len(labels)
                    regions.append(dict())
                labelIdx = labels[label]
                regions[labelIdx][name] = len(regions[labelIdx])


def sortMatrix(hm, regionsFileName, transcriptID, transcript_id_designator, verbose=True):
    """
    Iterate through the files noted by regionsFileName and sort hm accordingly
    """

    labels = dict()
    regions = []
    defaultGroup = None
    if len(regionsFileName) == 1:
        defaultGroup = "genes"
    for fname in regionsFileName:
        fp = dti.openPossiblyCompressed(fname)
        line = dti.getNext(fp)
        labelColumn = None
        while line.startswith("#"):
            if not labelColumn:
                labelColumn = dti.getLabel(line)
            line = dti.getNext(fp)
        while line.startswith("track "):
            line = dti.getNext(fp)

        # Find the label column
        subtract = 0
        if labelColumn is not None:
            subtract = 1

        # Determine the file type and load into a list (or list of lists)
        cols = line.strip().split("\t")
        if len(cols) - subtract < 3:
            raise RuntimeError('{0} does not seem to be a recognized file type!'.format(fname))
        elif len(cols) - subtract <= 6:
            loadBED(line, fp, fname, labelColumn, labels, regions, defaultGroup)
        elif len(cols) and dti.seemsLikeGTF(cols):
            loadGTF(line, fp, fname, labels, regions, transcriptID, transcript_id_designator, defaultGroup)
        else:
            loadBED(line, fp, fname, labelColumn, labels, regions, defaultGroup)
        fp.close()

    # Do some sanity checking on the group labels and region names within them
    s1 = set(hm.parameters['group_labels'])
    if verbose:
        for e in labels:
            if e not in s1:
                sys.exit("The computeMatrix output is missing the '{}' region group. It has {} but the specified regions have {}.\n".format(e, s1, labels.keys()))

    # Make a dictionary out of current labels and regions
    d = dict()
    pos = 0
    groupSizes = dict()
    for idx, label in enumerate(hm.parameters['group_labels']):
        s = hm.parameters['group_boundaries'][idx]
        e = hm.parameters['group_boundaries'][idx + 1]
        if label not in labels:
            continue
        d[label] = dict()
        groupSize = 0
        for reg in hm.matrix.regions[s:e]:
            d[label][reg[2]] = pos
            pos += 1
            groupSize += 1
        groupSizes[label] = groupSize

    # Convert labels to an ordered list
    labelsList = [""] * len(labels)
    for k, v in labels.items():
        labelsList[v] = k

    # Reorder
    order = []
    boundaries = [0]
    for idx, label in enumerate(labelsList):
        # Make an ordered list out of the region names in this region group
        _ = [""] * len(regions[idx])
        for k, v in regions[idx].items():
            _[v] = k
        sz = 0  # Track the number of enries actually matched
        for name in _:
            if name not in d[label]:
                if verbose:
                    sys.stderr.write("Skipping {}, due to being absent in the computeMatrix output.\n".format(name))
                continue
            sz += 1
            order.append(d[label][name])
        if sz == 0 and verbose:
            sys.exit("The region group {} had no matching entries!\n".format(label))
        boundaries.append(sz + boundaries[-1])
    hm.matrix.regions = [hm.matrix.regions[i] for i in order]
    order = np.array(order)
    hm.matrix.matrix = hm.matrix.matrix[order, :]

    # Update the parameters
    hm.parameters["group_labels"] = labelsList
    hm.matrix.group_labels = labelsList
    hm.parameters["group_boundaries"] = boundaries
    hm.matrix.group_boundaries = boundaries


def main(args=None):
    if len(sys.argv) == 1:
        args = ["-h"]
    if len(sys.argv) == 2:
        args = [sys.argv[1], "-h"]
    args = parse_arguments().parse_args(args)

    hm = heatmapper.heatmapper()
    if not isinstance(args.matrixFile, list):
        hm.read_matrix_file(args.matrixFile)
    if args.command == 'info':
        printInfo(hm)
    elif args.command == 'dataRange':
        printDataRange(hm)
    elif args.command == 'subset':
        sIdx = getSampleBounds(args, hm)
        gIdx, gBounds = getGroupBounds(args, hm)

        # groups
        hm.matrix.regions = subsetRegions(hm, gIdx)
        # matrix
        hm.matrix.matrix = hm.matrix.matrix[gIdx, :]
        hm.matrix.matrix = hm.matrix.matrix[:, sIdx]
        # boundaries
        if args.samples is None:
            args.samples = hm.matrix.sample_labels
        hm.matrix.sample_boundaries = hm.matrix.sample_boundaries[0:len(args.samples) + 1]
        hm.matrix.group_boundaries = gBounds.tolist()
        # special params
        keepIdx = set()
        for _, sample in enumerate(hm.matrix.sample_labels):
            if sample in args.samples:
                keepIdx.add(_)
        for param in hm.special_params:
            hm.parameters[param] = [v for k, v in enumerate(hm.parameters[param]) if k in keepIdx]
        # labels
        hm.matrix.sample_labels = args.samples
        if args.groups is None:
            args.groups = hm.matrix.group_labels
        hm.matrix.group_labels = args.groups
        # save
        hm.save_matrix(args.outFileName)
    elif args.command == 'filterStrand':
        filterHeatmap(hm, args)
        hm.save_matrix(args.outFileName)
    elif args.command == 'filterValues':
        filterHeatmapValues(hm, args.min, args.max)
        hm.save_matrix(args.outFileName)
    elif args.command == 'rbind':
        rbindMatrices(hm, args)
        hm.save_matrix(args.outFileName)
    elif args.command == 'cbind':
        cbindMatrices(hm, args)
        hm.save_matrix(args.outFileName)
    elif args.command == 'sort':
        sortMatrix(hm, args.regionsFileName, args.transcriptID, args.transcript_id_designator)
        hm.save_matrix(args.outFileName)
    elif args.command == 'relabel':
        relabelMatrix(hm, args)
        hm.save_matrix(args.outFileName)
    else:
        sys.exit("Unknown command {0}!\n".format(args.command))
