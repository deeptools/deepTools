import sys
import gzip
from collections import OrderedDict
import numpy as np
from copy import deepcopy

import pyBigWig
from deeptools import getScorePerBigWigBin
from deeptools import mapReduce
from deeptools.utilities import toString, toBytes, smartLabels
from deeptools.heatmapper_utilities import getProfileTicks


old_settings = np.seterr(all='ignore')


def chopRegions(exonsInput, left=0, right=0):
    """
    exons is a list of (start, end) tuples. The goal is to chop these into
    separate lists of tuples, to take care or unscaled regions. "left" and
    "right" denote regions of a given size to exclude from the normal binning
    process (unscaled regions).

    This outputs three lists of (start, end) tuples:

    leftBins: 5' unscaled regions
    bodyBins: body bins for scaling
    rightBins: 3' unscaled regions

    In addition are two integers
    padLeft: Number of bases of padding on the left (due to not being able to fulfill "left")
    padRight: As above, but on the right side
    """
    leftBins = []
    rightBins = []
    padLeft = 0
    padRight = 0
    exons = deepcopy(exonsInput)
    while len(exons) > 0 and left > 0:
        width = exons[0][1] - exons[0][0]
        if width <= left:
            leftBins.append(exons[0])
            del exons[0]
            left -= width
        else:
            leftBins.append((exons[0][0], exons[0][0] + left))
            exons[0] = (exons[0][0] + left, exons[0][1])
            left = 0
    if left > 0:
        padLeft = left

    while len(exons) > 0 and right > 0:
        width = exons[-1][1] - exons[-1][0]
        if width <= right:
            rightBins.append(exons[-1])
            del exons[-1]
            right -= width
        else:
            rightBins.append((exons[-1][1] - right, exons[-1][1]))
            exons[-1] = (exons[-1][0], exons[-1][1] - right)
            right = 0
    if right > 0:
        padRight = right

    return leftBins, exons, rightBins[::-1], padLeft, padRight


def chopRegionsFromMiddle(exonsInput, left=0, right=0):
    """
    Like chopRegions(), above, but returns two lists of tuples on each side of
    the center point of the exons.

    The steps are as follow:

     1) Find the center point of the set of exons (e.g., [(0, 200), (300, 400), (800, 900)] would be centered at 200)
       * If a given exon spans the center point then the exon is split
     2) The given number of bases at the end of the left-of-center list are extracted
       * If the set of exons don't contain enough bases, then padLeft is incremented accordingly
     3) As above but for the right-of-center list
     4) A tuple of (#2, #3, pading on the left, and padding on the right) is returned
    """
    leftBins = []
    rightBins = []
    size = sum([x[1] - x[0] for x in exonsInput])
    middle = size // 2
    cumulativeSum = 0
    padLeft = 0
    padRight = 0
    exons = deepcopy(exonsInput)

    # Split exons in half
    for exon in exons:
        size = exon[1] - exon[0]
        if cumulativeSum >= middle:
            rightBins.append(exon)
        elif cumulativeSum + size < middle:
            leftBins.append(exon)
        else:
            # Don't add 0-width exonic bins!
            if exon[0] < exon[1] - cumulativeSum - size + middle:
                leftBins.append((exon[0], exon[1] - cumulativeSum - size + middle))
            if exon[1] - cumulativeSum - size + middle < exon[1]:
                rightBins.append((exon[1] - cumulativeSum - size + middle, exon[1]))
        cumulativeSum += size

    # Trim leftBins/adjust padLeft
    lSum = sum([x[1] - x[0] for x in leftBins])
    if lSum > left:
        lSum = 0
        for i, exon in enumerate(leftBins[::-1]):
            size = exon[1] - exon[0]
            if lSum + size > left:
                leftBins[-i - 1] = (exon[1] + lSum - left, exon[1])
                break
            lSum += size
            if lSum == left:
                break
        i += 1
        if i < len(leftBins):
            leftBins = leftBins[-i:]
    elif lSum < left:
        padLeft = left - lSum

    # Trim rightBins/adjust padRight
    rSum = sum([x[1] - x[0] for x in rightBins])
    if rSum > right:
        rSum = 0
        for i, exon in enumerate(rightBins):
            size = exon[1] - exon[0]
            if rSum + size > right:
                rightBins[i] = (exon[0], exon[1] - rSum - size + right)
                break
            rSum += size
            if rSum == right:
                break
        rightBins = rightBins[:i + 1]
    elif rSum < right:
        padRight = right - rSum

    return leftBins, rightBins, padLeft, padRight


def trimZones(zones, maxLength, binSize, padRight):
    """
    Given a (variable length) list of lists of (start, end) tuples, trim/remove and tuple that extends past maxLength (e.g., the end of a chromosome)

    Returns the trimmed zones and padding
    """
    output = []
    for zone, nbins in zones:
        outZone = []
        changed = False
        for reg in zone:
            if reg[0] >= maxLength:
                changed = True
                padRight += reg[1] - reg[0]
                continue

            if reg[1] > maxLength:
                changed = True
                padRight += reg[1] - maxLength
                reg = (reg[0], maxLength)
            if reg[1] > reg[0]:
                outZone.append(reg)
        if changed:
            nBins = sum(x[1] - x[0] for x in outZone) // binSize
        else:
            nBins = nbins
        output.append((outZone, nBins))
    return output, padRight


def compute_sub_matrix_wrapper(args):
    return heatmapper.compute_sub_matrix_worker(*args)


class heatmapper(object):
    """
    Class to handle the reading and
    plotting of matrices.
    """

    def __init__(self):
        self.parameters = None
        self.lengthDict = None
        self.matrix = None
        self.regions = None
        self.blackList = None
        self.quiet = True
        # These are parameters that were single values in versions <3 but are now internally lists. See issue #614
        self.special_params = set(['unscaled 5 prime', 'unscaled 3 prime', 'body', 'downstream', 'upstream', 'ref point', 'bin size'])

    def getTicks(self, idx):
        """
        This is essentially a wrapper around getProfileTicks to accomdate the fact that each column has its own ticks.
        """
        xticks, xtickslabel = getProfileTicks(self, self.reference_point_label[idx], self.startLabel, self.endLabel, idx)
        return xticks, xtickslabel

    def computeMatrix(self, score_file_list, regions_file, parameters, blackListFileName=None, verbose=False, allArgs=None):
        """
        Splits into
        multiple cores the computation of the scores
        per bin for each region (defined by a hash '#'
        in the regions (BED/GFF) file.
        """
        if parameters['body'] > 0 and \
                parameters['body'] % parameters['bin size'] > 0:
            exit("The --regionBodyLength has to be "
                 "a multiple of --binSize.\nCurrently the "
                 "values are {} {} for\nregionsBodyLength and "
                 "binSize respectively\n".format(parameters['body'],
                                                 parameters['bin size']))

        # the beforeRegionStartLength is extended such that
        # length is a multiple of binSize
        if parameters['downstream'] % parameters['bin size'] > 0:
            exit("Length of region after the body has to be "
                 "a multiple of --binSize.\nCurrent value "
                 "is {}\n".format(parameters['downstream']))

        if parameters['upstream'] % parameters['bin size'] > 0:
            exit("Length of region before the body has to be a multiple of "
                 "--binSize\nCurrent value is {}\n".format(parameters['upstream']))

        if parameters['unscaled 5 prime'] % parameters['bin size'] > 0:
            exit("Length of the unscaled 5 prime region has to be a multiple of "
                 "--binSize\nCurrent value is {}\n".format(parameters['unscaled 5 prime']))

        if parameters['unscaled 3 prime'] % parameters['bin size'] > 0:
            exit("Length of the unscaled 5 prime region has to be a multiple of "
                 "--binSize\nCurrent value is {}\n".format(parameters['unscaled 3 prime']))

        if parameters['unscaled 5 prime'] + parameters['unscaled 3 prime'] > 0 and parameters['body'] == 0:
            exit('Unscaled 5- and 3-prime regions only make sense with the scale-regions subcommand.\n')

        # Take care of GTF options
        transcriptID = "transcript"
        exonID = "exon"
        transcript_id_designator = "transcript_id"
        keepExons = False
        self.quiet = False
        if allArgs is not None:
            allArgs = vars(allArgs)
            transcriptID = allArgs.get("transcriptID", transcriptID)
            exonID = allArgs.get("exonID", exonID)
            transcript_id_designator = allArgs.get("transcript_id_designator", transcript_id_designator)
            keepExons = allArgs.get("keepExons", keepExons)
            self.quiet = allArgs.get("quiet", self.quiet)

        chromSizes, _ = getScorePerBigWigBin.getChromSizes(score_file_list)
        res, labels = mapReduce.mapReduce([score_file_list, parameters],
                                          compute_sub_matrix_wrapper,
                                          chromSizes,
                                          self_=self,
                                          bedFile=regions_file,
                                          blackListFileName=blackListFileName,
                                          numberOfProcessors=parameters['proc number'],
                                          includeLabels=True,
                                          transcriptID=transcriptID,
                                          exonID=exonID,
                                          transcript_id_designator=transcript_id_designator,
                                          keepExons=keepExons,
                                          verbose=verbose)
        # each worker in the pool returns a tuple containing
        # the submatrix data, the regions that correspond to the
        # submatrix, and the number of regions lacking scores
        # Since this is largely unsorted, we need to sort by group

        # merge all the submatrices into matrix
        matrix = np.concatenate([r[0] for r in res], axis=0)
        regions = []
        regions_no_score = 0
        for idx in range(len(res)):
            if len(res[idx][1]):
                regions.extend(res[idx][1])
                regions_no_score += res[idx][2]
        groups = [x[3] for x in regions]
        foo = sorted(zip(groups, list(range(len(regions))), regions))
        sortIdx = [x[1] for x in foo]
        regions = [x[2] for x in foo]
        matrix = matrix[sortIdx]

        # mask invalid (nan) values
        matrix = np.ma.masked_invalid(matrix)

        assert matrix.shape[0] == len(regions), \
            "matrix length does not match regions length"

        if len(regions) == 0:
            sys.stderr.write("\nERROR: Either the BED file does not contain any valid regions or there are none remaining after filtering.\n")
            exit(1)
        if regions_no_score == len(regions):
            exit("\nERROR: None of the BED regions could be found in the bigWig"
                 "file.\nPlease check that the bigwig file is valid and "
                 "that the chromosome names between the BED file and "
                 "the bigWig file correspond to each other\n")

        if regions_no_score > len(regions) * 0.75:
            file_type = 'bigwig' if score_file_list[0].endswith(".bw") else "BAM"
            prcnt = 100 * float(regions_no_score) / len(regions)
            sys.stderr.write(
                "\n\nWarning: {0:.2f}% of regions are *not* associated\n"
                "to any score in the given {1} file. Check that the\n"
                "chromosome names from the BED file are consistent with\n"
                "the chromosome names in the given {2} file and that both\n"
                "files refer to the same species\n\n".format(prcnt,
                                                             file_type,
                                                             file_type))

        self.parameters = parameters

        numcols = matrix.shape[1]
        num_ind_cols = self.get_num_individual_matrix_cols()
        sample_boundaries = list(range(0, numcols + num_ind_cols, num_ind_cols))
        if allArgs is not None and allArgs['samplesLabel'] is not None:
            sample_labels = allArgs['samplesLabel']
        else:
            sample_labels = smartLabels(score_file_list)

        # Determine the group boundaries
        group_boundaries = []
        group_labels_filtered = []
        last_idx = -1
        for x in range(len(regions)):
            if regions[x][3] != last_idx:
                last_idx = regions[x][3]
                group_boundaries.append(x)
                group_labels_filtered.append(labels[last_idx])
        group_boundaries.append(len(regions))

        # check if a given group is too small. Groups that
        # are too small can't be plotted and an exception is thrown.
        group_len = np.diff(group_boundaries)
        if len(group_len) > 1:
            sum_len = sum(group_len)
            group_frac = [float(x) / sum_len for x in group_len]
            if min(group_frac) <= 0.002:
                sys.stderr.write(
                    "One of the groups defined in the bed file is "
                    "too small.\nGroups that are too small can't be plotted. "
                    "\n")

        self.matrix = _matrix(regions, matrix,
                              group_boundaries,
                              sample_boundaries,
                              group_labels_filtered,
                              sample_labels)

        if parameters['skip zeros']:
            self.matrix.removeempty()

    @staticmethod
    def compute_sub_matrix_worker(self, chrom, start, end, score_file_list, parameters, regions):
        """
        Returns
        -------
        numpy matrix
            A numpy matrix that contains per each row the values found per each of the regions given
        """
        if parameters['verbose']:
            sys.stderr.write("Processing {}:{}-{}\n".format(chrom, start, end))

        # read BAM or scores file
        score_file_handles = []
        for sc_file in score_file_list:
            score_file_handles.append(pyBigWig.open(sc_file))

        # determine the number of matrix columns based on the lengths
        # given by the user, times the number of score files
        matrix_cols = len(score_file_list) * \
            ((parameters['downstream'] +
              parameters['unscaled 5 prime'] + parameters['unscaled 3 prime'] +
              parameters['upstream'] + parameters['body']) //
             parameters['bin size'])

        # create an empty matrix to store the values
        sub_matrix = np.zeros((len(regions), matrix_cols))
        sub_matrix[:] = np.NAN

        j = 0
        sub_regions = []
        regions_no_score = 0
        for transcript in regions:
            feature_chrom = transcript[0]
            exons = transcript[1]
            feature_start = exons[0][0]
            feature_end = exons[-1][1]
            feature_name = transcript[2]
            feature_strand = transcript[4]
            padLeft = 0
            padRight = 0
            padLeftNaN = 0
            padRightNaN = 0
            upstream = []
            downstream = []

            # get the body length
            body_length = np.sum([x[1] - x[0] for x in exons]) - parameters['unscaled 5 prime'] - parameters['unscaled 3 prime']

            # print some information
            if parameters['body'] > 0 and \
                    body_length < parameters['bin size']:
                if not self.quiet:
                    sys.stderr.write("A region that is shorter than the bin size (possibly only after accounting for unscaled regions) was found: "
                                     "({0}) {1} {2}:{3}:{4}. Skipping...\n".format((body_length - parameters['unscaled 5 prime'] - parameters['unscaled 3 prime']),
                                                                                   feature_name, feature_chrom,
                                                                                   feature_start, feature_end))
                coverage = np.zeros(matrix_cols)
                if not parameters['missing data as zero']:
                    coverage[:] = np.nan
            else:
                if feature_strand == '-':
                    if parameters['downstream'] > 0:
                        upstream = [(feature_start - parameters['downstream'], feature_start)]
                    if parameters['upstream'] > 0:
                        downstream = [(feature_end, feature_end + parameters['upstream'])]
                    unscaled5prime, body, unscaled3prime, padLeft, padRight = chopRegions(exons, left=parameters['unscaled 3 prime'], right=parameters['unscaled 5 prime'])
                    # bins per zone
                    a = parameters['downstream'] // parameters['bin size']
                    b = parameters['unscaled 3 prime'] // parameters['bin size']
                    d = parameters['unscaled 5 prime'] // parameters['bin size']
                    e = parameters['upstream'] // parameters['bin size']
                else:
                    if parameters['upstream'] > 0:
                        upstream = [(feature_start - parameters['upstream'], feature_start)]
                    if parameters['downstream'] > 0:
                        downstream = [(feature_end, feature_end + parameters['downstream'])]
                    unscaled5prime, body, unscaled3prime, padLeft, padRight = chopRegions(exons, left=parameters['unscaled 5 prime'], right=parameters['unscaled 3 prime'])
                    a = parameters['upstream'] // parameters['bin size']
                    b = parameters['unscaled 5 prime'] // parameters['bin size']
                    d = parameters['unscaled 3 prime'] // parameters['bin size']
                    e = parameters['downstream'] // parameters['bin size']
                c = parameters['body'] // parameters['bin size']

                # build zones (each is a list of tuples)
                #  zone0: region before the region start,
                #  zone1: unscaled 5 prime region
                #  zone2: the body of the region
                #  zone3: unscaled 3 prime region
                #  zone4: the region from the end of the region downstream
                #  the format for each zone is: [(start, end), ...], number of bins
                # Note that for "reference-point", upstream/downstream will go
                # through the exons (if requested) and then possibly continue
                # on the other side (unless parameters['nan after end'] is true)
                if parameters['body'] > 0:
                    zones = [(upstream, a), (unscaled5prime, b), (body, c), (unscaled3prime, d), (downstream, e)]
                elif parameters['ref point'] == 'TES':  # around TES
                    if feature_strand == '-':
                        downstream, body, unscaled3prime, padRight, _ = chopRegions(exons, left=parameters['upstream'])
                        if padRight > 0 and parameters['nan after end'] is True:
                            padRightNaN += padRight
                        elif padRight > 0:
                            downstream.append((downstream[-1][1], downstream[-1][1] + padRight))
                        padRight = 0
                    else:
                        unscale5prime, body, upstream, _, padLeft = chopRegions(exons, right=parameters['upstream'])
                        if padLeft > 0 and parameters['nan after end'] is True:
                            padLeftNaN += padLeft
                        elif padLeft > 0:
                            upstream.insert(0, (upstream[0][0] - padLeft, upstream[0][0]))
                        padLeft = 0
                    e = np.sum([x[1] - x[0] for x in downstream]) // parameters['bin size']
                    a = np.sum([x[1] - x[0] for x in upstream]) // parameters['bin size']
                    zones = [(upstream, a), (downstream, e)]
                elif parameters['ref point'] == 'center':  # at the region center
                    if feature_strand == '-':
                        upstream, downstream, padLeft, padRight = chopRegionsFromMiddle(exons, left=parameters['downstream'], right=parameters['upstream'])
                    else:
                        upstream, downstream, padLeft, padRight = chopRegionsFromMiddle(exons, left=parameters['upstream'], right=parameters['downstream'])
                    if padLeft > 0 and parameters['nan after end'] is True:
                        padLeftNaN += padLeft
                    elif padLeft > 0:
                        if len(upstream) > 0:
                            upstream.insert(0, (upstream[0][0] - padLeft, upstream[0][0]))
                        else:
                            upstream = [(downstream[0][0] - padLeft, downstream[0][0])]
                    padLeft = 0
                    if padRight > 0 and parameters['nan after end'] is True:
                        padRightNaN += padRight
                    elif padRight > 0:
                        downstream.append((downstream[-1][1], downstream[-1][1] + padRight))
                    padRight = 0
                    a = np.sum([x[1] - x[0] for x in upstream]) // parameters['bin size']
                    e = np.sum([x[1] - x[0] for x in downstream]) // parameters['bin size']
                    # It's possible for a/e to be floats or 0 yet upstream/downstream isn't empty
                    if a < 1:
                        upstream = []
                        a = 0
                    if e < 1:
                        downstream = []
                        e = 0
                    zones = [(upstream, a), (downstream, e)]
                else:  # around TSS
                    if feature_strand == '-':
                        unscale5prime, body, upstream, _, padLeft = chopRegions(exons, right=parameters['downstream'])
                        if padLeft > 0 and parameters['nan after end'] is True:
                            padLeftNaN += padLeft
                        elif padLeft > 0:
                            upstream.insert(0, (upstream[0][0] - padLeft, upstream[0][0]))
                        padLeft = 0
                    else:
                        downstream, body, unscaled3prime, padRight, _ = chopRegions(exons, left=parameters['downstream'])
                        if padRight > 0 and parameters['nan after end'] is True:
                            padRightNaN += padRight
                        elif padRight > 0:
                            downstream.append((downstream[-1][1], downstream[-1][1] + padRight))
                        padRight = 0
                    a = np.sum([x[1] - x[0] for x in upstream]) // parameters['bin size']
                    e = np.sum([x[1] - x[0] for x in downstream]) // parameters['bin size']
                    zones = [(upstream, a), (downstream, e)]

                foo = parameters['upstream']
                bar = parameters['downstream']
                if feature_strand == '-':
                    foo, bar = bar, foo
                if padLeftNaN > 0:
                    expected = foo // parameters['bin size']
                    padLeftNaN = int(round(float(padLeftNaN) / parameters['bin size']))
                    if expected - padLeftNaN - a > 0:
                        padLeftNaN += 1
                if padRightNaN > 0:
                    expected = bar // parameters['bin size']
                    padRightNaN = int(round(float(padRightNaN) / parameters['bin size']))
                    if expected - padRightNaN - e > 0:
                        padRightNaN += 1

                coverage = []
                # compute the values for each of the files being processed.
                # "cov" is a numpy array of bins
                for sc_handler in score_file_handles:
                    # We're only supporting bigWig files at this point
                    cov = heatmapper.coverage_from_big_wig(
                        sc_handler, feature_chrom, zones,
                        parameters['bin size'],
                        parameters['bin avg type'],
                        parameters['missing data as zero'],
                        not self.quiet)

                    if padLeftNaN > 0:
                        cov = np.concatenate([[np.nan] * padLeftNaN, cov])
                    if padRightNaN > 0:
                        cov = np.concatenate([cov, [np.nan] * padRightNaN])

                    if feature_strand == "-":
                        cov = cov[::-1]

                    coverage = np.hstack([coverage, cov])

            if coverage is None:
                regions_no_score += 1
                if not self.quiet:
                    sys.stderr.write(
                        "No data was found for region "
                        "{0} {1}:{2}-{3}. Skipping...\n".format(
                            feature_name, feature_chrom,
                            feature_start, feature_end))

                coverage = np.zeros(matrix_cols)
                if not parameters['missing data as zero']:
                    coverage[:] = np.nan

            try:
                temp = coverage.copy()
                temp[np.isnan(temp)] = 0
            except:
                if not self.quiet:
                    sys.stderr.write(
                        "No scores defined for region "
                        "{0} {1}:{2}-{3}. Skipping...\n".format(feature_name,
                                                                feature_chrom,
                                                                feature_start,
                                                                feature_end))
                coverage = np.zeros(matrix_cols)
                if not parameters['missing data as zero']:
                    coverage[:] = np.nan

            if parameters['min threshold'] is not None and coverage.min() <= parameters['min threshold']:
                continue
            if parameters['max threshold'] is not None and coverage.max() >= parameters['max threshold']:
                continue
            if parameters['scale'] != 1:
                coverage = parameters['scale'] * coverage

            sub_matrix[j, :] = coverage

            sub_regions.append(transcript)
            j += 1

        # remove empty rows
        sub_matrix = sub_matrix[0:j, :]
        if len(sub_regions) != len(sub_matrix[:, 0]):
            sys.stderr.write("regions lengths do not match\n")
        return sub_matrix, sub_regions, regions_no_score

    @staticmethod
    def coverage_from_array(valuesArray, zones, binSize, avgType):
        try:
            valuesArray[0]
        except (IndexError, TypeError) as detail:
            sys.stderr.write("{0}\nvalues array value: {1}, zones {2}\n".format(detail, valuesArray, zones))

        cvglist = []
        zoneEnd = 0
        valStart = 0
        valEnd = 0
        for zone, nBins in zones:
            if nBins:
                # linspace is used to more or less evenly partition the data points into the given number of bins
                zoneEnd += nBins
                valStart = valEnd
                valEnd += np.sum([x[1] - x[0] for x in zone])
                counts_list = []

                # Partition the space into bins
                if nBins == 1:
                    pos_array = np.array([valStart])
                else:
                    pos_array = np.linspace(valStart, valEnd, nBins, endpoint=False, dtype=int)
                pos_array = np.append(pos_array, valEnd)

                idx = 0
                while idx < nBins:
                    idxStart = int(pos_array[idx])
                    idxEnd = max(int(pos_array[idx + 1]), idxStart + 1)
                    try:
                        counts_list.append(heatmapper.my_average(valuesArray[idxStart:idxEnd], avgType))
                    except Exception as detail:
                        sys.stderr.write("Exception found: {0}\n".format(detail))
                    idx += 1
                cvglist.append(np.array(counts_list))

        return np.concatenate(cvglist)

    @staticmethod
    def change_chrom_names(chrom):
        """
        Changes UCSC chromosome names to ensembl chromosome names
        and vice versa.
        """
        if chrom.startswith('chr'):
            # remove the chr part from chromosome name
            chrom = chrom[3:]
            if chrom == "M":
                chrom = "MT"
        else:
            # prefix with 'chr' the chromosome name
            chrom = 'chr' + chrom
            if chrom == "chrMT":
                chrom = "chrM"

        return chrom

    @staticmethod
    def coverage_from_big_wig(bigwig, chrom, zones, binSize, avgType, nansAsZeros=False, verbose=True):

        """
        uses pyBigWig
        to query a region define by chrom and zones.
        The output is an array that contains the bigwig
        value per base pair. The summary over bins is
        done in a later step when coverage_from_array is called.
        This method is more reliable than querying the bins
        directly from the bigwig, which should be more efficient.

        By default, any region, even if no chromosome match is found
        on the bigwig file, produces a result. In other words
        no regions are skipped.

        zones: array as follows zone0: region before the region start,
                                zone1: 5' unscaled region (if present)
                                zone2: the body of the region (not always present)
                                zone3: 3' unscaled region (if present)
                                zone4: the region from the end of the region downstream

               each zone is a tuple containing start, end, and number of bins


        This is useful if several matrices wants to be merged
        or if the sorted BED output of one computeMatrix operation
        needs to be used for other cases
        """
        nVals = 0
        for zone, _ in zones:
            for region in zone:
                nVals += region[1] - region[0]

        values_array = np.zeros(nVals)
        if not nansAsZeros:
            values_array[:] = np.nan
        if chrom not in list(bigwig.chroms().keys()):
            unmod_name = chrom
            chrom = heatmapper.change_chrom_names(chrom)
            if chrom not in list(bigwig.chroms().keys()):
                if verbose:
                    sys.stderr.write("Warning: Your chromosome names do not match.\nPlease check that the "
                                     "chromosome names in your BED file\ncorrespond to the names in your "
                                     "bigWig file.\nAn empty line will be added to your heatmap.\nThe problematic "
                                     "chromosome name is {0}\n\n".format(unmod_name))

                # return empty nan array
                return heatmapper.coverage_from_array(values_array, zones, binSize, avgType)

        maxLen = bigwig.chroms(chrom)
        startIdx = 0
        endIdx = 0
        for zone, _ in zones:
            for region in zone:
                startIdx = endIdx
                if region[0] < 0:
                    endIdx += abs(region[0])
                    values_array[startIdx:endIdx] = np.nan
                    startIdx = endIdx
                start = max(0, region[0])
                end = min(maxLen, region[1])
                endIdx += end - start
                if start < end:
                    # This won't be the case if we extend off the front of a chromosome, such as (-100, 0)
                    values_array[startIdx:endIdx] = bigwig.values(chrom, start, end)
                if end < region[1]:
                    startIdx = endIdx
                    endIdx += region[1] - end
                    values_array[startIdx:endIdx] = np.nan

        # replaces nans for zeros
        if nansAsZeros:
            values_array[np.isnan(values_array)] = 0

        return heatmapper.coverage_from_array(values_array, zones,
                                              binSize, avgType)

    @staticmethod
    def my_average(valuesArray, avgType='mean'):
        """
        computes the mean, median, etc but only for those values
        that are not Nan
        """
        valuesArray = np.ma.masked_invalid(valuesArray)
        avg = np.ma.__getattribute__(avgType)(valuesArray)
        if isinstance(avg, np.ma.core.MaskedConstant):
            return np.nan
        else:
            return avg

    def matrix_from_dict(self, matrixDict, regionsDict, parameters):
        self.regionsDict = regionsDict
        self.matrixDict = matrixDict
        self.parameters = parameters
        self.lengthDict = OrderedDict()
        self.matrixAvgsDict = OrderedDict()

    def read_matrix_file(self, matrix_file):
        # reads a bed file containing the position
        # of genomic intervals
        # In case a hash sign '#' is found in the
        # file, this is considered as a delimiter
        # to split the heatmap into groups

        import json
        regions = []
        matrix_rows = []
        current_group_index = 0
        max_group_bound = None

        fh = gzip.open(matrix_file)
        for line in fh:
            line = toString(line).strip()
            # read the header file containing the parameters
            # used
            if line.startswith("@"):
                # the parameters used are saved using
                # json
                self.parameters = json.loads(line[1:].strip())
                max_group_bound = self.parameters['group_boundaries'][1]
                continue

            # split the line into bed interval and matrix values
            region = line.split('\t')
            chrom, start, end, name, score, strand = region[0:6]
            matrix_row = np.ma.masked_invalid(np.fromiter(region[6:], np.float64))
            matrix_rows.append(matrix_row)
            starts = start.split(",")
            ends = end.split(",")
            regs = [(int(x), int(y)) for x, y in zip(starts, ends)]
            # get the group index
            if len(regions) >= max_group_bound:
                current_group_index += 1
                max_group_bound = self.parameters['group_boundaries'][current_group_index + 1]
            regions.append([chrom, regs, name, max_group_bound, strand, score])

        matrix = np.vstack(matrix_rows)
        self.matrix = _matrix(regions, matrix, self.parameters['group_boundaries'],
                              self.parameters['sample_boundaries'],
                              group_labels=self.parameters['group_labels'],
                              sample_labels=self.parameters['sample_labels'])

        if 'sort regions' in self.parameters:
            self.matrix.set_sorting_method(self.parameters['sort regions'],
                                           self.parameters['sort using'])

        # Versions of computeMatrix before 3.0 didn't have an entry of these per column, fix that
        nSamples = len(self.matrix.sample_labels)
        h = dict()
        for k, v in self.parameters.items():
            if k in self.special_params and type(v) is not list:
                v = [v] * nSamples
                if len(v) == 0:
                    v = [None] * nSamples
            h[k] = v
        self.parameters = h

        return

    def save_matrix(self, file_name):
        """
        saves the data required to reconstruct the matrix
        the format is:
        A header containing the parameters used to create the matrix
        encoded as:
        @key:value\tkey2:value2 etc...
        The rest of the file has the same first 5 columns of a
        BED file: chromosome name, start, end, name, score and strand,
        all separated by tabs. After the fifth column the matrix
        values are appended separated by tabs.
        Groups are separated by adding a line starting with a hash (#)
        and followed by the group name.

        The file is gzipped.
        """
        import json
        self.parameters['sample_labels'] = self.matrix.sample_labels
        self.parameters['group_labels'] = self.matrix.group_labels
        self.parameters['sample_boundaries'] = self.matrix.sample_boundaries
        self.parameters['group_boundaries'] = self.matrix.group_boundaries

        # Redo the parameters, ensuring things related to ticks and labels are repeated appropriately
        nSamples = len(self.matrix.sample_labels)
        h = dict()
        for k, v in self.parameters.items():
            if type(v) is list and len(v) == 0:
                v = None
            if k in self.special_params and type(v) is not list:
                v = [v] * nSamples
                if len(v) == 0:
                    v = [None] * nSamples
            h[k] = v
        fh = gzip.open(file_name, 'wb')
        params_str = json.dumps(h, separators=(',', ':'))
        fh.write(toBytes("@" + params_str + "\n"))
        score_list = np.ma.masked_invalid(np.mean(self.matrix.matrix, axis=1))
        for idx, region in enumerate(self.matrix.regions):
            # join np_array values
            # keeping nans while converting them to strings
            if not np.ma.is_masked(score_list[idx]):
                np.float64(score_list[idx])
            matrix_values = "\t".join(
                np.char.mod('%f', self.matrix.matrix[idx, :]))
            starts = ["{0}".format(x[0]) for x in region[1]]
            ends = ["{0}".format(x[1]) for x in region[1]]
            starts = ",".join(starts)
            ends = ",".join(ends)
            # BEDish format (we don't currently store the score)
            fh.write(
                toBytes('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(
                        region[0],
                        starts,
                        ends,
                        region[2],
                        region[5],
                        region[4],
                        matrix_values)))
        fh.close()

    def save_tabulated_values(self, file_handle, reference_point_label='TSS', start_label='TSS', end_label='TES', averagetype='mean'):
        """
        Saves the values averaged by col using the avg_type
        given

        Args:
            file_handle: file name to save the file
            reference_point_label: Name of the reference point label
            start_label: Name of the star label
            end_label: Name of the end label
            averagetype: average type (e.g. mean, median, std)

        """
        #  get X labels
        w = self.parameters['bin size']
        b = self.parameters['upstream']
        a = self.parameters['downstream']
        c = self.parameters.get('unscaled 5 prime', 0)
        d = self.parameters.get('unscaled 3 prime', 0)
        m = self.parameters['body']

        xticks = []
        xtickslabel = []
        for idx in range(self.matrix.get_num_samples()):
            if b[idx] < 1e5:
                quotient = 1000
                symbol = 'Kb'
            else:
                quotient = 1e6
                symbol = 'Mb'

            if m[idx] == 0:
                last = 0
                if len(xticks):
                    last = xticks[-1]
                xticks.extend([last + (k / w[idx]) for k in [w[idx], b[idx], b[idx] + a[idx]]])
                xtickslabel.extend(['{0:.1f}{1}'.format(-(float(b[idx]) / quotient), symbol), reference_point_label,
                                    '{0:.1f}{1}'.format(float(a[idx]) / quotient, symbol)])

            else:
                xticks_values = [w[idx]]

                # only if upstream region is set, add a x tick
                if b[idx] > 0:
                    xticks_values.append(b[idx])
                    xtickslabel.append('{0:.1f}{1}'.format(-(float(b[idx]) / quotient), symbol))

                xtickslabel.append(start_label)

                if c[idx] > 0:
                    xticks_values.append(b[idx] + c[idx])
                    xtickslabel.append("")

                if d[idx] > 0:
                    xticks_values.append(b[idx] + c[idx] + m[idx])
                    xtickslabel.append("")

                xticks_values.append(b[idx] + c[idx] + m[idx] + d[idx])
                xtickslabel.append(end_label)

                if a[idx] > 0:
                    xticks_values.append(b[idx] + c[idx] + m[idx] + d[idx] + a[idx])
                    xtickslabel.append('{0:.1f}{1}'.format(float(a[idx]) / quotient, symbol))

                last = 0
                if len(xticks):
                    last = xticks[-1]
                xticks.extend([last + (k / w[idx]) for k in xticks_values])
        x_axis = np.arange(xticks[-1]) + 1
        labs = []
        for x_value in x_axis:
            if x_value in xticks and xtickslabel[xticks.index(x_value)]:
                labs.append(xtickslabel[xticks.index(x_value)])
            elif x_value in xticks:
                labs.append("tick")
            else:
                labs.append("")

        with open(file_handle, 'w') as fh:
            # write labels
            fh.write("bin labels\t\t{}\n".format("\t".join(labs)))
            fh.write('bins\t\t{}\n'.format("\t".join([str(x) for x in x_axis])))

            for sample_idx in range(self.matrix.get_num_samples()):
                for group_idx in range(self.matrix.get_num_groups()):
                    sub_matrix = self.matrix.get_matrix(group_idx, sample_idx)
                    values = [str(x) for x in np.ma.__getattribute__(averagetype)(sub_matrix['matrix'], axis=0)]
                    fh.write("{}\t{}\t{}\n".format(sub_matrix['sample'], sub_matrix['group'], "\t".join(values)))

    def save_matrix_values(self, file_name):
        # print a header telling the group names and their length
        fh = open(file_name, 'wb')
        info = []
        groups_len = np.diff(self.matrix.group_boundaries)
        for i in range(len(self.matrix.group_labels)):
            info.append("{}:{}".format(self.matrix.group_labels[i],
                                       groups_len[i]))
        fh.write(toBytes("#{}\n".format("\t".join(info))))
        # add to header the x axis values
        fh.write(toBytes("#downstream:{}\tupstream:{}\tbody:{}\tbin size:{}\tunscaled 5 prime:{}\tunscaled 3 prime:{}\n".format(
                 self.parameters['downstream'],
                 self.parameters['upstream'],
                 self.parameters['body'],
                 self.parameters['bin size'],
                 self.parameters.get('unscaled 5 prime', 0),
                 self.parameters.get('unscaled 3 prime', 0))))
        sample_len = np.diff(self.matrix.sample_boundaries)
        for i in range(len(self.matrix.sample_labels)):
            info.extend([self.matrix.sample_labels[i]] * sample_len[i])
        fh.write(toBytes("{}\n".format("\t".join(info))))

        fh.close()
        # reopen again using append mode
        fh = open(file_name, 'ab')
        np.savetxt(fh, self.matrix.matrix, fmt="%.4g", delimiter="\t")
        fh.close()

    def save_BED(self, file_handle):
        boundaries = np.array(self.matrix.group_boundaries)
        # Add a header
        file_handle.write("#chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRGB\tblockCount\tblockSizes\tblockStart\tdeepTools_group")
        if self.matrix.silhouette is not None:
            file_handle.write("\tsilhouette")
        file_handle.write("\n")
        for idx, region in enumerate(self.matrix.regions):
            # the label id corresponds to the last boundary
            # that is smaller than the region index.
            # for example for a boundary array = [0, 10, 20]
            # and labels ['a', 'b', 'c'],
            # for index 5, the label is 'a', for
            # index 10, the label is 'b' etc
            label_idx = np.flatnonzero(boundaries <= idx)[-1]
            starts = ["{0}".format(x[0]) for x in region[1]]
            ends = ["{0}".format(x[1]) for x in region[1]]
            starts = ",".join(starts)
            ends = ",".join(ends)
            file_handle.write(
                '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{1}\t{2}\t0'.format(
                    region[0],
                    region[1][0][0],
                    region[1][-1][1],
                    region[2],
                    region[5],
                    region[4]))
            file_handle.write(
                '\t{0}\t{1}\t{2}\t{3}'.format(
                    len(region[1]),
                    ",".join([str(int(y) - int(x)) for x, y in region[1]]),
                    ",".join([str(int(x) - int(starts[0])) for x, y in region[1]]),
                    self.matrix.group_labels[label_idx]))
            if self.matrix.silhouette is not None:
                file_handle.write("\t{}".format(self.matrix.silhouette[idx]))
            file_handle.write("\n")
        file_handle.close()

    @staticmethod
    def matrix_avg(matrix, avgType='mean'):
        matrix = np.ma.masked_invalid(matrix)
        return np.ma.__getattribute__(avgType)(matrix, axis=0)

    def get_individual_matrices(self, matrix):
        """In case multiple matrices are saved one after the other
        this method splits them appart.
        Returns a list containing the matrices
        """
        num_cols = matrix.shape[1]
        num_ind_cols = self.get_num_individual_matrix_cols()
        matrices_list = []
        for i in range(0, num_cols, num_ind_cols):
            if i + num_ind_cols > num_cols:
                break
            matrices_list.append(matrix[:, i:i + num_ind_cols])
        return matrices_list

    def get_num_individual_matrix_cols(self):
        """
        returns the number of columns  that
        each matrix should have. This is done because
        the final matrix that is plotted can be composed
        of smaller matrices that are merged one after
        the other.
        """
        matrixCols = ((self.parameters['downstream'] + self.parameters['upstream'] + self.parameters['body'] + self.parameters['unscaled 5 prime'] + self.parameters['unscaled 3 prime']) //
                      self.parameters['bin size'])

        return matrixCols


def computeSilhouetteScore(d, idx, labels):
    """
    Given a square distance matrix with NaN diagonals, compute the silhouette score
    of a given row (idx). Each row should have an associated label (labels).
    """
    keep = ~np.isnan(d[idx, ])
    foo = np.bincount(labels[keep], weights=d[idx, ][keep])
    groupSizes = np.bincount(labels[keep])
    intraIdx = labels[idx]
    if groupSizes[intraIdx] == 1:
        return 0
    intra = foo[labels[idx]] / groupSizes[intraIdx]
    interMask = np.arange(len(foo))[np.arange(len(foo)) != labels[idx]]
    inter = np.min(foo[interMask] / groupSizes[interMask])
    return (inter - intra) / max(inter, intra)


class _matrix(object):
    """
    class to hold heatmapper matrices
    The base data is a large matrix
    with definition to know the boundaries for row and col divisions.
    Col divisions represent groups within a subset, e.g. Active and
    inactive from PolII bigwig data.

    Row division represent different samples, for example
    PolII in males vs. PolII in females.

    This is an internal class of the heatmapper class
    """

    def __init__(self, regions, matrix, group_boundaries, sample_boundaries,
                 group_labels=None, sample_labels=None):

        # simple checks
        assert matrix.shape[0] == group_boundaries[-1], \
            "row max do not match matrix shape"
        assert matrix.shape[1] == sample_boundaries[-1], \
            "col max do not match matrix shape"

        self.regions = regions
        self.matrix = matrix
        self.group_boundaries = group_boundaries
        self.sample_boundaries = sample_boundaries
        self.sort_method = None
        self.sort_using = None
        self.silhouette = None

        if group_labels is None:
            self.group_labels = ['group {}'.format(x)
                                 for x in range(len(group_boundaries) - 1)]
        else:
            assert len(group_labels) == len(group_boundaries) - 1, \
                "number of group labels does not match number of groups"
            self.group_labels = group_labels

        if sample_labels is None:
            self.sample_labels = ['sample {}'.format(x)
                                  for x in range(len(sample_boundaries) - 1)]
        else:
            assert len(sample_labels) == len(sample_boundaries) - 1, \
                "number of sample labels does not match number of samples"
            self.sample_labels = sample_labels

    def get_matrix(self, group, sample):
        """
        Returns a sub matrix from the large
        matrix. Group and sample are ids,
        thus, row = 0, col=0 get the first group
        of the first sample.

        Returns
        -------
        dictionary containing the matrix,
        the group label and the sample label
        """
        group_start = self.group_boundaries[group]
        group_end = self.group_boundaries[group + 1]
        sample_start = self.sample_boundaries[sample]
        sample_end = self.sample_boundaries[sample + 1]

        return {'matrix': np.ma.masked_invalid(self.matrix[group_start:group_end, :][:, sample_start:sample_end]),
                'group': self.group_labels[group],
                'sample': self.sample_labels[sample]}

    def get_num_samples(self):
        return len(self.sample_labels)

    def get_num_groups(self):
        return len(self.group_labels)

    def set_group_labels(self, new_labels):
        """ sets new labels for groups
        """
        if len(new_labels) != len(self.group_labels):
            raise ValueError("length new labels != length original labels")
        self.group_labels = new_labels

    def set_sample_labels(self, new_labels):
        """ sets new labels for groups
        """
        if len(new_labels) != len(self.sample_labels):
            raise ValueError("length new labels != length original labels")
        self.sample_labels = new_labels

    def set_sorting_method(self, sort_method, sort_using):
        self.sort_method = sort_method
        self.sort_using = sort_using

    def get_regions(self):
        """Returns the regions per group

        Returns
        ------
        list

            Each element of the list is itself a list
            of dictionaries containing the regions info:
            chrom, start, end, strand, name etc.

            Each element of the list corresponds to each
            of the groups
        """
        regions = []
        for idx in range(len(self.group_labels)):
            start = self.group_boundaries[idx]
            end = self.group_boundaries[idx + 1]
            regions.append(self.regions[start:end])

        return regions

    def sort_groups(self, sort_using='mean', sort_method='no', sample_list=None):
        """
        Sorts and rearranges the submatrices according to the
        sorting method given.
        """
        if sort_method == 'no':
            return

        if (sample_list is not None) and (len(sample_list) > 0):
            # get the ids that correspond to the selected sample list
            idx_to_keep = []
            for sample_idx in sample_list:
                idx_to_keep += range(self.sample_boundaries[sample_idx], self.sample_boundaries[sample_idx + 1])

            matrix = self.matrix[:, idx_to_keep]

        else:
            matrix = self.matrix

        # compute the row average:
        if sort_using == 'region_length':
            matrix_avgs = list()
            for x in self.regions:
                matrix_avgs.append(np.sum([bar[1] - bar[0] for bar in x[1]]))
            matrix_avgs = np.array(matrix_avgs)
        elif sort_using == 'mean':
            matrix_avgs = np.nanmean(matrix, axis=1)
        elif sort_using == 'mean':
            matrix_avgs = np.nanmean(matrix, axis=1)
        elif sort_using == 'median':
            matrix_avgs = np.nanmedian(matrix, axis=1)
        elif sort_using == 'max':
            matrix_avgs = np.nanmax(matrix, axis=1)
        elif sort_using == 'min':
            matrix_avgs = np.nanmin(matrix, axis=1)
        elif sort_using == 'sum':
            matrix_avgs = np.nansum(matrix, axis=1)
        else:
            sys.exit("{} is an unsupported sorting method".format(sort_using))

        # order per group
        _sorted_regions = []
        _sorted_matrix = []
        for idx in range(len(self.group_labels)):
            start = self.group_boundaries[idx]
            end = self.group_boundaries[idx + 1]
            order = matrix_avgs[start:end].argsort()
            if sort_method == 'descend':
                order = order[::-1]
            _sorted_matrix.append(self.matrix[start:end, :][order, :])
            # sort the regions
            _reg = self.regions[start:end]
            for idx in order:
                _sorted_regions.append(_reg[idx])

        self.matrix = np.vstack(_sorted_matrix)
        self.regions = _sorted_regions
        self.set_sorting_method(sort_method, sort_using)

    def hmcluster(self, k, evaluate_silhouette=True, method='kmeans', clustering_samples=None):
        matrix = np.asarray(self.matrix)
        matrix_to_cluster = matrix
        if clustering_samples is not None:
            assert all(i > 0 for i in clustering_samples),\
                "all indices should be bigger than or equal to 1."
            assert all(i <= len(self.sample_labels) for i in
                       clustering_samples),\
                "each index should be smaller than or equal to {}(total "\
                "number of samples.)".format(len(self.sample_labels))

            clustering_samples = np.asarray(clustering_samples) - 1

            samples_cols = []
            for idx in clustering_samples:
                samples_cols += range(self.sample_boundaries[idx],
                                      self.sample_boundaries[idx + 1])

            matrix_to_cluster = matrix_to_cluster[:, samples_cols]
        if np.any(np.isnan(matrix_to_cluster)):
            # replace nans for 0 otherwise kmeans produces a weird behaviour
            sys.stderr.write("*Warning* For clustering nan values have to be replaced by zeros \n")
            matrix_to_cluster[np.isnan(matrix_to_cluster)] = 0

        if method == 'kmeans':
            from scipy.cluster.vq import vq, kmeans

            centroids, _ = kmeans(matrix_to_cluster, k)
            # order the centroids in an attempt to
            # get the same cluster order
            cluster_labels, _ = vq(matrix_to_cluster, centroids)

        if method == 'hierarchical':
            # normally too slow for large data sets
            from scipy.cluster.hierarchy import fcluster, linkage
            Z = linkage(matrix_to_cluster, method='ward', metric='euclidean')
            cluster_labels = fcluster(Z, k, criterion='maxclust')
            # hierarchical clustering labels from 1 .. k
            # while k-means labels 0 .. k -1
            # Thus, for consistency, we subtract 1
            cluster_labels -= 1

        # sort clusters
        _clustered_mean = []
        _cluster_ids_list = []
        for cluster in range(k):
            cluster_ids = np.flatnonzero(cluster_labels == cluster)
            _cluster_ids_list.append(cluster_ids)
            _clustered_mean.append(matrix_to_cluster[cluster_ids, :].mean())

        # reorder clusters based on mean
        cluster_order = np.argsort(_clustered_mean)[::-1]
        # create groups using the clustering
        self.group_labels = []
        self.group_boundaries = [0]
        _clustered_regions = []
        _clustered_matrix = []
        cluster_number = 1
        for cluster in cluster_order:
            self.group_labels.append("cluster_{}".format(cluster_number))
            cluster_number += 1
            cluster_ids = _cluster_ids_list[cluster]
            self.group_boundaries.append(self.group_boundaries[-1] +
                                         len(cluster_ids))
            _clustered_matrix.append(self.matrix[cluster_ids, :])
            for idx in cluster_ids:
                _clustered_regions.append(self.regions[idx])

        self.regions = _clustered_regions
        self.matrix = np.vstack(_clustered_matrix)

        return idx

    def computeSilhouette(self, k):
        if k > 1:
            from scipy.spatial.distance import pdist, squareform

            silhouette = np.repeat(0.0, self.group_boundaries[-1])
            groupSizes = np.subtract(self.group_boundaries[1:], self.group_boundaries[:-1])
            labels = np.repeat(np.arange(k), groupSizes)

            d = pdist(self.matrix)
            d2 = squareform(d)
            np.fill_diagonal(d2, np.nan)  # This excludes the diagonal
            for idx in range(len(labels)):
                silhouette[idx] = computeSilhouetteScore(d2, idx, labels)
            sys.stderr.write("The average silhouette score is: {}\n".format(np.mean(silhouette)))
            self.silhouette = silhouette

    def removeempty(self):
        """
        removes matrix rows containing only zeros or nans
        """
        to_keep = []
        score_list = np.ma.masked_invalid(np.mean(self.matrix, axis=1))
        for idx, region in enumerate(self.regions):
            if np.ma.is_masked(score_list[idx]) or np.float64(score_list[idx]) == 0:
                continue
            else:
                to_keep.append(idx)
        self.regions = [self.regions[x] for x in to_keep]
        self.matrix = self.matrix[to_keep, :]
        # adjust sample boundaries
        to_keep = np.array(to_keep)
        self.group_boundaries = [len(to_keep[to_keep < x]) for x in self.group_boundaries]

    def flatten(self):
        """
        flatten and remove nans from matrix. Useful
        to get max and mins from matrix.

        :return flattened matrix
        """
        matrix_flatten = np.asarray(self.matrix.flatten())
        # nans are removed from the flattened array
        matrix_flatten = matrix_flatten[~np.isnan(matrix_flatten)]
        if len(matrix_flatten) == 0:
            num_nan = len(np.flatnonzero(np.isnan(self.matrix.flatten())))
            raise ValueError("matrix only contains nans "
                             "(total nans: {})".format(num_nan))
        return matrix_flatten
