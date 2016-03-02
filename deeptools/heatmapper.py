import sys
from os.path import splitext, basename
import gzip
from collections import OrderedDict
import numpy as np
import multiprocessing

# NGS packages
import pysam

import pyBigWig
import deeptools.readBed
from deeptools import mapReduce


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

    def computeMatrix(self, score_file_list, regions_file, parameters, blackListFileName=None, verbose=False):
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

        regions, group_labels = self.get_regions_and_groups(regions_file, blackListFileName=blackListFileName, verbose=verbose)

        # args to pass to the multiprocessing workers
        mp_args = []
        # prepare groups of regions to send to workers.
        regions_per_worker = 400 / len(score_file_list)
        for index in range(0, len(regions), regions_per_worker):
            index_end = min(len(regions), index + regions_per_worker)
            mp_args.append((score_file_list, regions[index:index_end],
                            parameters))

        if len(mp_args) > 1 and parameters['proc number'] > 1:
            pool = multiprocessing.Pool(parameters['proc number'])
            res = pool.map_async(compute_sub_matrix_wrapper,
                                 mp_args).get(9999999)
        else:
            res = map(compute_sub_matrix_wrapper, mp_args)

        # each worker in the pool returns a tuple containing
        # the submatrix data and the regions that correspond to the
        # submatrix

        # merge all the submatrices into matrix
        matrix = np.concatenate([r[0] for r in res], axis=0)
        regions = []
        regions_no_score = 0
        for idx in range(len(res)):
            if len(res[idx][1]):
                regions.extend(res[idx][1])
                regions_no_score += res[idx][2]

        # mask invalid (nan) values
        matrix = np.ma.masked_invalid(matrix)

        assert matrix.shape[0] == len(regions), \
            "matrix length does not match regions length"

        if len(regions) == 0:
            sys.stderr.write(
                "\nERROR: BED file does not contain any valid regions. "
                "Please check\n")
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
                "\n\nWarning: {:.2f}% of regions are *not* associated\n"
                "to any score in the given {} file. Check that the\n"
                "chromosome names from the BED file are consistent with\n"
                "the chromosome names in the given {} file and that both\n"
                "files refer to the same species\n\n".format(prcnt,
                                                             file_type,
                                                             file_type))

        self.parameters = parameters

        numcols = matrix.shape[1]
        num_ind_cols = self.get_num_individual_matrix_cols()
        sample_boundaries = range(0, numcols + num_ind_cols, num_ind_cols)
        sample_labels = [splitext(basename(x))[0] for x in score_file_list]

        # Determine the group boundaries, since any filtering out of regions will change things
        group_boundaries = []
        group_labels_filtered = []
        last_idx = -1
        for x in range(len(regions)):
            if regions[x].group_idx != last_idx:
                last_idx = regions[x].group_idx
                group_boundaries.append(x)
                group_labels_filtered.append(group_labels[last_idx])
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
    def compute_sub_matrix_worker(score_file_list, regions,
                                  parameters):

        """
        Parameters
        ----------
        score_file_list : list
            List of strings. Contains the list of
            bam or bigwig files to be used.
        regions : list of dictionaries
            Each item in the list is a dictionary
            containing containing the fields: 'chrom', 'start', 'end', 'name' and 'strand.
        parameters : dict
            Contains values that specify the length of bins, the number of bp after a reference point etc.

        Returns
        -------
        numpy matrix
            A numpy matrix that contains per each row the values found per each of the regions given
        """

        # read BAM or scores file
        score_file_handlers = []
        for sc_file in score_file_list:
            if sc_file.endswith(".bam"):
                score_file_handlers.append(pysam.Samfile(sc_file, 'rb'))
            else:
                score_file_handlers.append(pyBigWig.open(sc_file))

        """
        if score_file_list[0].endswith(".bam"):
            bamfile_list = []
            for score_file in score_file_list:
                bamfile_list.append(pysam.Samfile(score_file, 'rb'))
        else:
            bigwig_list = []
            from bx.bbi.bigwig_file import BigWigFile
            for score_file in score_file_list:
                bigwig_list.append(BigWigFile(file=open(score_file, 'r' )))
        """

        # determine the number of matrix columns based on the lengths
        # given by the user, times the number of score files
        matrix_cols = len(score_file_list) * \
            ((parameters['downstream'] +
              parameters['unscaled 5 prime'] + parameters['unscaled 3 prime'] +
              parameters['upstream'] + parameters['body']) /
             parameters['bin size'])

        # create an empty matrix to store the values
        sub_matrix = np.zeros((len(regions), matrix_cols))
        sub_matrix[:] = np.NAN

        j = 0
        sub_regions = []
        regions_no_score = 0
        for feature in regions:
            # print some information
            if parameters['body'] > 0 and \
                    feature.end - feature.start - parameters['unscaled 5 prime'] - parameters['unscaled 3 prime'] < parameters['bin size']:
                if parameters['verbose']:
                    sys.stderr.write("A region that is shorter than the bin size (possibly only after accounting for unscaled regions) was found: "
                                     "({}) {} {}:{}:{}. Skipping...\n".format((feature.end - feature.start - parameters['unscaled 5 prime'] - parameters['unscaled 3 prime']),
                                                                              feature.name, feature.chrom,
                                                                              feature.start, feature.end))
                coverage = np.zeros(matrix_cols)
                coverage[:] = np.nan
            else:
                if feature.strand == '-':
                    a = parameters['upstream'] / parameters['bin size']
                    b = parameters['downstream'] / parameters['bin size']
                    d = parameters['unscaled 5 prime'] / parameters['bin size']
                    c = parameters['unscaled 3 prime'] / parameters['bin size']
                    start = feature.end - parameters['unscaled 5 prime']
                    end = feature['start'] + parameters['unscaled 3 prime']
                else:
                    b = parameters['upstream'] / parameters['bin size']
                    a = parameters['downstream'] / parameters['bin size']
                    c = parameters['unscaled 5 prime'] / parameters['bin size']
                    d = parameters['unscaled 3 prime'] / parameters['bin size']
                    start = feature['start'] + parameters['unscaled 5 prime']
                    end = feature.end - parameters['unscaled 3 prime']

                # build zones:
                #  zone0: region before the region start,
                #  zone1: unscaled 5 prime region
                #  zone2: the body of the region (not always present)
                #  zone3: unscaled 3 prime region
                #  zone4: the region from the end of the region downstream
                #  the format for each zone is: start, end, number of bins
                if parameters['body'] > 0:
                    zones = \
                        [(feature.start - b * parameters['bin size'], feature.start, b),
                         (feature.start, feature.start + c * parameters['bin size'], c),
                         (feature.start + c * parameters['bin size'], feature.end - d * parameters['bin size'], parameters['body'] / parameters['bin size']),
                         (feature.end - d * parameters['bin size'], feature.end, d),
                         (feature.end, feature.end + a * parameters['bin size'], a)]
                elif parameters['ref point'] == 'TES':  # around TES
                    zones = [(end - b * parameters['bin size'], end, b),
                             (end, end + a * parameters['bin size'], a)]
                elif parameters['ref point'] == 'center':  # at the region center
                    middle_point = feature.start + (feature.end - feature.start) / 2
                    zones = [(middle_point - b * parameters['bin size'],
                              middle_point, b),
                             (middle_point,
                              middle_point + a * parameters['bin size'], a)]
                else:  # around TSS
                    zones = [(start - b * parameters['bin size'], start, b),
                             (start, start + a * parameters['bin size'], a)]

                if feature.start - b * parameters['bin size'] < 0:
                    if parameters['verbose']:
                        sys.stderr.write(
                            "Warning:region too close to chromosome start "
                            "for {} {}:{}:{}.\n".format(feature.name,
                                                        feature.chrom,
                                                        feature.start,
                                                        feature.end))
                # check if after extension, the region extends beyond the
                # chromosome length
                if feature.chrom not in score_file_handlers[0].chroms().keys():
                    chrom = heatmapper.change_chrom_names(feature.chrom)
                else:
                    chrom = feature.chrom
                if feature.end + a * parameters['bin size'] > score_file_handlers[0].chroms(chrom):
                    if parameters['verbose']:
                        sys.stderr.write(
                            "Warning:region extends beyond chromosome end "
                            "for {} {}:{}:{}.\n".format(feature.name,
                                                        feature.chrom,
                                                        feature.start,
                                                        feature.end))

                coverage = []
                # compute the values (coverage in the case of bam files)
                # for each of the files being processed.
                for idx, sc_handler in enumerate(score_file_handlers):
                    # check if the file is bam or bigwig
                    if score_file_list[idx].endswith(".bam"):
                        cov = heatmapper.coverage_from_bam(
                            sc_handler, feature.chrom, zones,
                            parameters['bin size'],
                            parameters['bin avg type'],
                            parameters['verbose'])
                    else:
                        cov = heatmapper.coverage_from_big_wig(
                            sc_handler, feature.chrom, zones,
                            parameters['bin size'],
                            parameters['bin avg type'],
                            parameters['missing data as zero'],
                            parameters['verbose'])

                    if feature.strand == "-":
                        cov = cov[::-1]

                    if parameters['nan after end'] and parameters['body'] == 0 \
                            and parameters['ref point'] == 'TSS':
                        # convert the gene length to bin length
                        region_length_in_bins = (feature.end - feature.start) / parameters['bin size']
                        b = parameters['upstream'] / parameters['bin size']
                        # convert to nan any region after the end of the region
                        cov[b + region_length_in_bins:] = np.nan

                    coverage = np.hstack([coverage, cov])

                """
                else:
                    for bigwig in bigwig_list:
                        cov = heatmapper.coverage_from_big_wig(
                            bigwig, feature.chrom, zones,
                            parameters['bin size'],
                            parameters['bin avg type'],
                            parameters['missing data as zero'])
                        if feature.strand == "-":
                            cov = cov[::-1]
                        coverage = np.hstack([coverage, cov])
                """

            if coverage is None:
                regions_no_score += 1
                if parameters['verbose']:
                    sys.stderr.write(
                        "No data was found for region "
                        "{} {}:{}-{}. Skipping...\n".format(
                            feature.name, feature.chrom,
                            feature.start, feature.end))

                coverage = np.zeros(matrix_cols)
                if not parameters['missing data as zero']:
                    coverage[:] = np.nan

            try:
                temp = coverage.copy()
                temp[np.isnan(temp)] = 0
            except:
                if parameters['verbose']:
                    sys.stderr.write(
                        "No scores defined for region "
                        "{} {}:{}-{}. Skipping...\n".format(feature.name,
                                                            feature.chrom,
                                                            feature.start,
                                                            feature.end))
                coverage = np.zeros(matrix_cols)
                if not parameters['missing data as zero']:
                    coverage[:] = np.nan

            if parameters['min threshold'] and coverage.min() <= parameters['min threshold']:
                continue
            if parameters['max threshold'] and coverage.max() >= parameters['max threshold']:
                continue
            if parameters['scale'] != 1:
                coverage = parameters['scale'] * coverage

            sub_matrix[j, :] = coverage

            sub_regions.append(feature)
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
            sys.stderr.write("{}\nvalues array value: {}, zones {}\n".format(detail, valuesArray, zones))

        cvglist = []
        start = zones[0][0]
        for zone_start, zone_end, num_bins in zones:
                # the linspace is to get equally spaced positions along the range
            # If the gene is short the sampling regions could overlap,
            # if it is long, the sampling regions would be spaced
            counts_list = []

            # this case happens when the downstream or upstream
            # region is set to 0
            if zone_start == zone_end:
                continue
            if num_bins == 1:
                pos_array, step_size = np.array([zone_start]), zone_end - zone_start
            else:
                (pos_array, step_size) = np.linspace(zone_start, zone_end, num_bins,
                                                     endpoint=False,
                                                     retstep=True)
            step_size = np.ceil(step_size)

            for pos in np.floor(pos_array):
                index_start = int(pos - start)
                index_end = int(index_start + step_size)
                try:
                    counts_list.append(
                        heatmapper.my_average(valuesArray[index_start:index_end],
                                              avgType))
                except Exception as detail:
                    sys.stderr.write("Exception found. "
                                     "Message: {}\n".format(detail))
            cvglist.append(np.array(counts_list))
        return np.concatenate(cvglist)

    @staticmethod
    def change_chrom_names(chrom):
        """
        Changes UCSC chromosome names to ensembl chromosome names
        and vice versa.
        """
        # TODO: mapping from chromosome names, e.g., mt, unknown
        if chrom.startswith('chr'):
            # remove the chr part from chromosome name
            chrom = chrom[3:]
        else:
            # prefix with 'chr' the chromosome name
            chrom = 'chr' + chrom

        return chrom

    @staticmethod
    def coverage_from_bam(bamfile, chrom, zones, binSize, avgType, verbose=True):
        """
        currently this method is deactivated because is too slow.
        It is preferred to create a coverage bigiwig file from the
        bam file and then run heatmapper.
        """
        if chrom not in bamfile.references:
            chrom = heatmapper.change_chrom_names(chrom)
            if chrom not in bamfile.references:
                sys.stderr.write(
                    "Skipping region located at unknown chromosome: {} "
                    "Known chromosomes are: {}\n".format(chrom,
                                                         bamfile.references))
                return None
            elif verbose:
                sys.stderr.write("Warning: Your chromosome names do "
                                 "not match.\n Please check that the "
                                 "chromosome names in your BED "
                                 "file correspond to the names in your "
                                 "bigWig file.\n An empty line will be "
                                 "added to your heatmap."
                                 "scheme.\n")

        start = zones[0][0]
        end = zones[-1][1]
        try:
            values_array = np.zeros(end - start)
            for read in bamfile.fetch(chrom, min(0, start), end):
                index_start = max(read.pos - start, 0)
                index_end = min(read.pos - start + read.qlen, end - start)
                values_array[index_start:index_end] += 1
        except ValueError:
            sys.stderr.write("Value out of range for region {}s {} {}\n".format(chrom, start, end))
            return np.array([0])  # return something innocuous

        return heatmapper.coverage_from_array(values_array, zones, binSize, avgType)

    @staticmethod
    def coverage_from_big_wig(bigwig, chrom, zones, binSize, avgType, nansAsZeros=False, verbose=True):

        """
        uses bigwig file reader from bx-python
        to query a region define by chrom and zones.
        The output is an array that contains the bigwig
        value per base pair. The summary over bins is
        done in a later step when coverage_from_array is called.
        This method is more reliable than quering the bins
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

        # initialize values array. The length of the array
        # is the length of the region which is defined
        # by the start of the first zone zones[0][0]
        # to the end of the last zone zones[-1][1]
        values_array = np.zeros(zones[-1][1] - zones[0][0])
        if not nansAsZeros:
            values_array[:] = np.nan
        bw_array = None
        if chrom not in bigwig.chroms().keys():
            unmod_name = chrom
            chrom = heatmapper.change_chrom_names(chrom)
            if chrom not in bigwig.chroms().keys():
                if verbose:
                    sys.stderr.write("Warning: Your chromosome names do not match.\nPlease check that the "
                                     "chromosome names in your BED file\ncorrespond to the names in your "
                                     "bigWig file.\nAn empty line will be added to your heatmap.\nThe problematic "
                                     "chromosome name is {}\n\n".format(unmod_name))

                # return empty nan array
                return heatmapper.coverage_from_array(values_array, zones, binSize, avgType)
        try:

            start = max(0, zones[0][0])
            end = min(bigwig.chroms(chrom), zones[-1][1])
            bw_array = bigwig.values(chrom, start, end)
        except Exception as detail:
                sys.stderr.write("Exception found. Message: "
                                 "{}\n".format(detail))
                sys.stderr.write("Problematic region: {}:{}-{}\n".format(chrom, zones[-1][1], zones[0][0]))

        # adjust bw_array if it extends beyond chromosome limits
        if bw_array is not None:
            if zones[0][0] < 0 and zones[-1][1] > bigwig.chroms(chrom):
                values_array[abs(zones[0][0]):len(bw_array) + abs(zones[0][0])] = bw_array
            elif zones[0][0] < 0:
                values_array[abs(zones[0][0]):] = bw_array
            elif zones[-1][1] > bigwig.chroms(chrom):
                values_array[:len(bw_array)] = bw_array
            else:
                values_array = np.array(bw_array)

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
        avg = np.__getattribute__(avgType)(valuesArray)
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

        fh = gzip.open(matrix_file)
        for line in fh:
            line = line.strip()
            # read the header file containing the parameters
            # used
            if line.startswith("@"):
                # the parameters used are saved using
                # json
                self.parameters = json.loads(line[1:].strip())
                continue

            # split the line into bed interval and matrix values
            region = line.split('\t')
            chrom, start, end, name, score, strand = region[0:6]
            matrix_row = np.ma.masked_invalid(np.fromiter(region[6:], np.float))
            matrix_rows.append(matrix_row)
            regions.append({'chrom': chrom, 'start': int(start),
                            'end': int(end), 'name': name, 'score': score,
                            'strand': strand})

        matrix = np.vstack(matrix_rows)
        self.matrix = _matrix(regions, matrix, self.parameters['group_boundaries'],
                              self.parameters['sample_boundaries'],
                              group_labels=self.parameters['group_labels'],
                              sample_labels=self.parameters['sample_labels'])

        if 'sort regions' in self.parameters:
            self.matrix.set_sorting_method(self.parameters['sort regions'],
                                           self.parameters['sort using'])
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

        fh = gzip.open(file_name, 'wb')
        params_str = json.dumps(self.parameters, separators=(',', ':'))
        fh.write("@" + params_str + "\n")
        score_list = np.ma.masked_invalid(np.mean(self.matrix.matrix, axis=1))
        for idx, region in enumerate(self.matrix.regions):
            # join np_array values
            # keeping nans while converting them to strings
            if not np.ma.is_masked(score_list[idx]):
                np.float(score_list[idx])
            matrix_values = "\t".join(
                np.char.mod('%f', self.matrix.matrix[idx, :]))
            fh.write(
                '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    region['chrom'],
                    region['start'],
                    region['end'],
                    region['name'],
                    region['score'],
                    region['strand'],
                    matrix_values))
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
        c = self.parameters['unscaled 5 prime']
        d = self.parameters['unscaled 3 prime']
        m = self.parameters['body']

        if b < 1e5:
            quotient = 1000
            symbol = 'Kb'
        else:
            quotient = 1e6
            symbol = 'Mb'

        if m == 0:
            xticks = [(k / w) for k in [w, b, b + a]]
            xtickslabel = ['{0:.1f}{1}'.format(-(float(b) / quotient), symbol), reference_point_label,
                           '{0:.1f}{1}'.format(float(a) / quotient, symbol)]

        else:
            xticks_values = [w]
            xtickslabel = []

            # only if upstream region is set, add a x tick
            if b > 0:
                xticks_values.append(b)
                xtickslabel.append('{0:.1f}{1}'.format(-(float(b) / quotient), symbol))

            xtickslabel.append(start_label)

            if c > 0:
                xticks_values.append(b + c)
                xtickslabel.append("")

            if d > 0:
                xticks_values.append(b + c + m)
                xtickslabel.append("")

            xticks_values.append(b + c + m + d)
            xtickslabel.append(end_label)

            if a > 0:
                xticks_values.append(b + c + m + d + a)
                xtickslabel.append('{0:.1f}{1}'.format(float(a) / quotient, symbol))

            xticks = [(k / w) for k in xticks_values]
        x_axis = np.arange(xticks[-1]) + 1
        labs = []
        for x_value in x_axis:
            if x_value in xticks:
                labs.append(xtickslabel[xticks.index(x_value)])
            else:
                labs.append("")

        with open(file_handle, 'w') as fh:
            # write labels
            fh.write("bin labels\t\t{}\n".format("\t".join(labs)))
            fh.write('bins\t\t{}\n'.format("\t".join([str(x) for x in x_axis])))

            for sample_idx in range(self.matrix.get_num_samples()):
                for group_idx in range(self.matrix.get_num_groups()):
                    sub_matrix = self.matrix.get_matrix(group_idx, sample_idx)
                    values = [str(x) for x in np.__getattribute__(averagetype)(sub_matrix['matrix'], axis=0)]
                    fh.write("{}\t{}\t{}\n".format(sub_matrix['sample'], sub_matrix['group'], "\t".join(values)))

    def save_matrix_values(self, file_name):
        # print a header telling the group names and their length
        fh = open(file_name, 'w')
        info = []
        groups_len = np.diff(self.matrix.group_boundaries)
        for i in range(len(self.matrix.group_labels)):
            info.append("{}:{}".format(self.matrix.group_labels[i],
                                       groups_len[i]))
        fh.write("#{}\n".format("\t".join(info)))
        # add to header the x axis values
        fh.write("#downstream:{}\tupstream:{}\tbody:{}\tbin size:{}\tunscaled 5 prime:{}\tunscaled 3 prime:{}\n".format(
                 self.parameters['downstream'],
                 self.parameters['upstream'],
                 self.parameters['body'],
                 self.parameters['bin size'],
                 self.parameters['unscaled 5 prime'],
                 self.parameters['unscaled 3 prime']))

        fh.close()
        # reopen again using append mode
        fh = open(file_name, 'a')
        np.savetxt(fh, self.matrix.matrix, fmt="%.4g")
        fh.close()

    def save_BED(self, file_handle):
        boundaries = np.array(self.matrix.group_boundaries)
        for idx, region in enumerate(self.matrix.regions):
            # the label id corresponds to the last boundary
            # that is smaller than the region index.
            # for example for a boundary array = [0, 10, 20]
            # and labels ['a', 'b', 'c'],
            # for index 5, the label is 'a', for
            # index 10, the label is 'b' etc
            label_idx = np.flatnonzero(boundaries <= idx)[-1]
            file_handle.write(
                '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    region['chrom'],
                    region['start'],
                    region['end'],
                    region['name'],
                    region['score'],
                    region['strand'],
                    self.matrix.group_labels[label_idx]))
            if idx + 1 in boundaries:
                file_handle.write('#{}\n'.format(
                    self.matrix.group_labels[label_idx]))
        file_handle.close()

    @staticmethod
    def matrix_avg(matrix, avgType='mean'):
        matrix = np.ma.masked_invalid(matrix)
        return np.__getattribute__(avgType)(matrix, axis=0)

    @staticmethod
    def get_regions_and_groups(regions_file, onlyMultiplesOf=1,
                               default_group_name='genes',
                               blackListFileName=None,
                               verbose=None):
        """
        Reads a bed file.
        In case is hash sign '#' is found in the
        file, this is considered as a delimiter
        to split the heatmap into groups

        Returns a list of regions with a label
        index appended to each and a list of labels
        """

        regions = []
        previnterval = None
        duplicates = 0
        totalintervals = 0
        groupintervals = 0
        includedintervals = 0
        group_labels = []
        group_idx = 0
        bed_file = deeptools.readBed.ReadBed(regions_file)
        blackList = None
        if blackListFileName is not None:
            blackList = mapReduce.BED_to_interval_tree(open(blackListFileName, "r"))

        for ginterval in bed_file:
            if ginterval.line.startswith("track") or ginterval.line.startswith("browser"):
                continue

            if ginterval.line.startswith('#'):
                # check for labels with no associated entries
                if groupintervals == 0:
                    continue
                else:
                    groupintervals = 0
                group_idx += 1
                label = ginterval.line[1:].strip()
                if label in group_labels:
                    # loop to find a unique label name
                    i = 0
                    while True:
                        i += 1
                        newlabel = label + "_r" + str(i)
                        if newlabel not in group_labels:
                            break

                group_labels.append(label)
                continue

            # Exclude blacklist overlaps
            if mapReduce.blOverlap(blackList, ginterval.chrom, [ginterval.start, ginterval.end]):
                continue

            # if the list of regions is to big, only
            # consider a fraction of the data
            # if totalintervals % onlyMultiplesOf != 0:
            #    continue
            # check for regions that have the same position as the previous.
            # This assumes that the regions file given is sorted
            totalintervals += 1
            if previnterval is not None:
                if previnterval.chrom == ginterval.chrom and \
                   previnterval.start == ginterval.start and \
                   previnterval.end == ginterval.end and \
                   previnterval.strand == ginterval.strand:
                    if verbose:
                        try:
                            genename = ginterval.name
                        except:
                            genename = ''
                        sys.stderr.write("*Warning* Duplicated region: "
                                         "{} {}:{}-{}.\n".format(genename,
                                                                 ginterval.chrom,
                                                                 ginterval.start,
                                                                 ginterval.end))
                    duplicates += 1

            groupintervals += 1
            previnterval = ginterval
            ginterval.group_idx = group_idx
            regions.append(ginterval)
            includedintervals += 1

        # in case we reach the end of the file
        # without encountering a hash,
        # a default name is given to regions
        using_default_group_name = False
        if not group_labels:
            group_labels.append(default_group_name)
            using_default_group_name = True
            groupintervals = 0

        # If there are any remaining intervals with no group label then add a fake one
        if groupintervals > 0:
            # There was a missing "#" at the end
            label = default_group_name
            if label in group_labels:
                # loop to find a unique label name
                i = 0
                while True:
                    i += 1
                    newlabel = label + "_r" + str(i)
                    if newlabel not in group_labels:
                        break
            group_labels.append(label)

        if verbose and duplicates > 0:
            sys.stderr.write(
                "{} ({:.2f}) regions covering the exact same interval "
                "were found".format(duplicates,
                                    float(duplicates) * 100 / totalintervals))

        if verbose:
            sys.stderr.write("Found:\n\tintervals: {}\n".format(len(regions)))
            if using_default_group_name:
                sys.stderr.write("\tno groups found\n")
            else:
                sys.stderr.write("\tgroups: {} [{}]\n\n".format(len(group_labels), ", ".join(group_labels)))

        return regions, group_labels

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
        matrixCols = ((self.parameters['downstream'] + self.parameters['upstream'] + self.parameters['body'] + self.parameters['unscaled 5 prime'] + self.parameters['unscaled 3 prime']) /
                      self.parameters['bin size'])

        return matrixCols


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

    def sort_groups(self, sort_using='mean', sort_method='no'):
        """
        Sorts and rearranges the submatrices according to the
        sorting method given.
        """
        if sort_method == 'no':
            return

        # compute the row average:
        if sort_using == 'region_length':
            matrix_avgs = np.array([x['end'] - x['start']
                                   for x in self.regions])
        else:
            matrix_avgs = np.__getattribute__(sort_using)(
                self.matrix, axis=1)

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

    def hmcluster(self, k, method='kmeans'):

        matrix = np.asarray(self.matrix)
        if np.any(np.isnan(matrix)):
            # replace nans for 0 otherwise kmeans produces a weird behaviour
            sys.stderr.write("*Warning* For clustering nan values have to be replaced by zeros \n")
            matrix[np.isnan(matrix)] = 0

        if method == 'kmeans':
            from scipy.cluster.vq import vq, kmeans

            centroids, _ = kmeans(matrix, k)
            # order the centroids in an attempt to
            # get the same cluster order
            order = np.argsort(centroids.mean(axis=1))[::-1]
            cluster_labels, _ = vq(matrix, centroids[order, :])

        if method == 'hierarchical':
            # normally too slow for large data sets
            from scipy.cluster.hierarchy import fcluster, linkage
            Z = linkage(matrix, method='ward', metric='euclidean')
            cluster_labels = fcluster(Z, k, criterion='maxclust')

        # create groups using the clustering
        self.group_labels = []
        self.group_boundaries = [0]
        _clustered_regions = []
        _clustered_matrix = []
        for cluster in range(k):
            self.group_labels.append("cluster {}".format(cluster + 1))
            cluster_ids = np.flatnonzero(cluster_labels == cluster)
            self.group_boundaries.append(self.group_boundaries[-1] +
                                         len(cluster_ids))
            _clustered_matrix.append(self.matrix[cluster_ids, :])
            for idx in cluster_ids:
                _clustered_regions.append(self.regions[idx])

        self.regions = _clustered_regions
        self.matrix = np.vstack(_clustered_matrix)
        return idx

    def removeempty(self):
        """
        removes matrix rows containing only zeros or nans
        """
        to_keep = []
        score_list = np.ma.masked_invalid(np.mean(self.matrix, axis=1))
        for idx, region in enumerate(self.regions):
            if np.ma.is_masked(score_list[idx]) or np.float(score_list[idx]) == 0:
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
