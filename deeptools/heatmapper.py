import sys
import gzip
from collections import OrderedDict
import numpy as np
import multiprocessing

# NGS packages
import pysam
from bx.intervals.io import GenomicIntervalReader


def compute_sub_matrix_wrapper(args):
    return heatmapper.compute_sub_matrix_worker(*args)


class heatmapper:
    """
    Class to handle the reading and
    plotting of matrices.
    """

    def __init__(self):
        self.regionsDict = None
        self.matrixDict = None
        self.parameters = None
        self.lengthDict = None
        self.matrixAvgsDict = None

    def computeMatrix(self, score_file, regions_file, parameters,
                      verbose=False):
        """
        Splits into
        multiple cores the computation of the scores
        per bin for each region (defined by a hash '#'
        in the regions (BED/GFF) file.
        """
        if parameters['body'] > 0 and \
                parameters['body'] % parameters['bin size'] > 0:
            sys.stderr.write("The --regionBodyLength has to be "
                             "a multiple of --binSize.\nCurrently the "
                             "values are {} {} for\nregionsBodyLength and "
                             "binSize respectively\n".format(
                    parameters['body'],
                    parameters['bin size']))
            exit(1)

        # the beforeRegionStartLength is extended such that
        # length is a multiple of binSize
        if parameters['downstream'] % parameters['bin size'] > 0:
            sys.stderr.write(
                "Length of region after the body has to be "
                "a multiple of --binSize.\nCurrent value "
                "is {}\n".format(parameters['downstream']))
            exit(1)

        if parameters['upstream'] % parameters['bin size'] > 0:
            sys.stderr.write(
                "Length of region before the body has to be a multiple of "
                "--binSize\nCurrent value is {}\n".format(
                    parameters['upstream']))
            exit(1)

        # determine the number of matrix columns based on the lengths
        # given by the user
        matrixCols = ((parameters['downstream'] +
                       parameters['upstream'] + parameters['body']) /
                      parameters['bin size'])

        regionsDict = self.getRegionsAndGroups(regions_file, verbose=verbose)
        group_len = [len(x) for x in regionsDict]
        # check if a given group is too small. Groups that
        # are too small can be plotted and an error is shown.
        if len(group_len) > 1:
            sum_len = sum(group_len)
            group_frac = [float(x)/sum_len for x in group_len]
            if min(group_frac) <= 0.002:
                sys.stderr.write(
                    "One of the groups defined in the bed file is "
                    "too small.\nGroups that are too small can't be plotted. "
                    "Please remove the group to continue.\n")
                exit(1)
        matrixDict = OrderedDict()
        matrixAvgsDict = OrderedDict()
        for label, regions in regionsDict.iteritems():
            # args to pass to the multiprocessing workers
            mp_args = []

            # prepare groups of 400 regions to send to workers.
            for index in range(0, len(regions), 400):
                index_end = min(len(regions), index + 400 )
                mp_args.append((score_file, regions[index:index_end],
                                matrixCols, parameters))

            if len(mp_args) > 1 and parameters['proc number'] > 1:
                pool = multiprocessing.Pool(parameters['proc number'])
                res = pool.map_async(compute_sub_matrix_wrapper,
                                     mp_args).get(9999999)
            else:
                res = map(compute_sub_matrix_wrapper, mp_args)

            # each worker in the pools returns a tuple containing
            # the submatrix data and the regions that correspond to the
            # submatrix

            # merge all the submatrices into matrix
            matrix = np.concatenate([r[0] for r in res], axis=0)
            # mask invalid (nan) values
            matrix = np.ma.masked_invalid(matrix)
            # merge all valid regions
            regionList = np.concatenate([r[1] for r in res], axis=0)

            regions_no_score = sum([r[2] for r in res])
            if len(regions) == 0:
                sys.stderr.write(
                    "\nERROR: BED file does not contain any valid regions. "
                    "Please check\n")
                exit(1)
            if regions_no_score == len(regions):
                sys.stderr.write(
                    "\nERROR: None of the BED regions could be found in the bigWig"
                     "file.\nPlease check that the bigwig file is valid and "
                     "that the chromosome names between the BED file and "
                     "the bigWig file correspond to each other\n")
                exit(1)
            if regions_no_score > len(regions) * 0.75 or len(regionList) == 0:
                file_type = 'bigwig' if score_file.endswith(".bw") else "BAM"
                prcnt = 100 * float(regions_no_score) / len(regions)
                sys.stderr.write(
                    "\n\nWarning: {:.2f}% of regions are *not* associated\n"
                    "to any score in the given {} file. Check that the\n"
                    "chromosome names from the BED file are consistent with\n"
                    "the chromosome names in the given {} file and that both\n"
                    "files refer to the same species\n\n".format(prcnt,
                                                                 file_type,
                                                                 file_type))

            matrixAvgsDict[label] = np.mean(matrix, axis=1)
            matrixDict[label] = matrix
            regionsDict[label] = regionList

        self.matrixDict = matrixDict
        self.regionsDict = regionsDict
        self.parameters = parameters
        self.matrixAvgsDict = matrixAvgsDict

    @staticmethod
    def compute_sub_matrix_worker(score_file, regions, matrixCols, parameters):
        # read BAM or scores file
        if score_file.endswith(".bam"):
            bamfile = pysam.Samfile(score_file, 'rb')
        else:
            from bx.bbi.bigwig_file import BigWigFile
            bigwig = BigWigFile(file=open(score_file, 'r' ))
        # create an empty matrix to store the values
        subMatrix = np.zeros((len(regions), matrixCols))
        subMatrix[:] = np.NAN

        j = 0
        subRegions = []
        regions_no_score = 0
        for feature in regions:
           # print some information
            if parameters['body'] > 0 and \
                    feature['end'] - feature['start'] < parameters['bin size']:
                if parameters['verbose']:
                    sys.stderr.write("A region that is shorter than "
                                     "then bin size was found: "
                                     "({}) {} {}:{}:{}. Skipping...\n".format(
                            (feature['end'] - feature['start']),
                            feature['name'], feature['chrom'],
                            feature['start'], feature['end']))
                continue

            if feature['strand'] == '-':
                a = parameters['upstream'] / parameters['bin size']
                b = parameters['downstream']  / parameters['bin size']
                start = feature['end']
                end = feature['start']
            else:
                b = parameters['upstream'] / parameters['bin size']
                a = parameters['downstream'] / parameters['bin size']
                start = feature['start']
                end = feature['end']

            # build zones:
            #  zone0: region before the region start,
            #  zone1: the body of the region (not always present)
            #  zone2: the region from the end of the region downstream
            #  the format for each zone is: start, end, number of bins
            if parameters['body'] > 0:
                zones = \
                    [(feature['start'] - b * parameters['bin size'],
                      feature['start'], b ),
                     (feature['start'],
                      feature['end'],
                      #feature['end'] - parameters['body'] /
                      #parameters['bin size'],
                      parameters['body'] / parameters['bin size']),
                     (feature['end'],
                      feature['end'] + a * parameters['bin size'], a)]
            elif parameters['ref point'] == 'TES':  # around TES
                zones = [(end - b * parameters['bin size'], end, b ),
                         (end, end + a * parameters['bin size'], a )]
            elif parameters['ref point'] == 'center':  # at the region center
                middlePoint = feature['start'] + (feature['end'] -
                                                  feature['start']) / 2
                zones = [(middlePoint - b * parameters['bin size'],
                          middlePoint, b),
                         (middlePoint,
                          middlePoint + a * parameters['bin size'], a)]
            else:  # around TSS
                zones = [(start - b * parameters['bin size'], start, b ),
                         (start, start + a * parameters['bin size'], a )]

            if feature['start'] - b * parameters['bin size'] < 0:
                if parameters['verbose']:
                    sys.stderr.write(
                        "Warning:region too close to chromosome start "
                        "for {} {}:{}:{}.\n".format(feature['name'],
                                                   feature['chrom'],
                                                   feature['start'],
                                                   feature['end']))
            coverage = None
            if score_file.endswith(".bam"):
                coverage = heatmapper.coverageFromBam(
                    bamfile, feature['chrom'], zones,
                    parameters['bin size'],
                    parameters['bin avg type'])

            else:
                coverage = heatmapper.coverageFromBigWig(
                    bigwig, feature['chrom'], zones,
                    parameters['bin size'],
                    parameters['bin avg type'],
                    parameters['missing data as zero'])

            """ 
            if coverage is None:
                regions_no_score += 1
                if parameters['verbose']:
                    sys.stderr.write(
                        "No data was found for region "
                        "{} {}:{}-{}. Skipping...\n".format(
                            feature['name'], feature['chrom'],
                            feature['start'], feature['end']))

                coverage = np.zeros(matrixCols)
                if not parameters['missing data as zero']:
                    coverage[:] = np.nan
                continue
            """
            try:
                temp = coverage.copy()
                temp[np.isnan(temp)] = 0
                totalScore = np.sum(temp)
            except:
                if parameters['verbose']:
                    sys.stderr.write(
                        "No scores defined for region "
                        "{} {}:{}-{}. Skipping...\n".format(feature['name'],
                                                            feature['chrom'],
                                                            feature['start'],
                                                            feature['end']))
                coverage = np.zeros(matrixCols)
                if not parameters['missing data as zero']:
                    coverage[:] = np.nan
                # to induce skipping if zero regions are omited this
                # variable is set to zero
                totalScore = 0

            if totalScore == 0:
                regions_no_score += 1
                if parameters['skip zeros']:
                    if parameters['verbose']:
                        sys.stderr.write(
                            "Skipping region with all scores equal to zero "
                            "for\n'{}' {}:{}-{}.\n\n".format(feature['name'],
                                                             feature['chrom'],
                                                             feature['start'],
                                                             feature['end']))
                    continue
                elif parameters['verbose']:
                    sys.stderr.write(
                        "Warning: All values are zero for "
                        "{} {}:{}-{}.\n".format(feature['name'],
                                                feature['chrom'],
                                                feature['start'],
                                                feature['end']))
                    sys.stderr.write(
                        "add --skipZeros to exclude such regions\n")

            if parameters['min threshold'] and \
                    coverage.min() <= parameters['min threshold']:
                continue
            if parameters['max threshold'] and \
                    coverage.max() >= parameters['max threshold']:
                continue
            if parameters['scale'] != 1:
                coverage = parameters['scale'] * coverage

            if feature['strand'] == "-":
                subMatrix[j, :] = coverage[::-1]
            else:
                subMatrix[j, :] = coverage

            if parameters['nan after end'] and parameters['body'] == 0 \
                    and parameters['ref point'] == 'TSS':
                # convert the gene length to bin length
                region_length_in_bins = \
                    (feature['end'] - feature['start']) / \
                    parameters['bin size']
                b = parameters['upstream'] / parameters['bin size']
                # convert to nan any region after the end of the region
                subMatrix[j, b + region_length_in_bins:] = np.nan

            subRegions.append(feature)
            j += 1

        # remove empty rows
        subMatrix = subMatrix[0:j, :]
        if len(subRegions) != len(subMatrix[:, 0]):
            sys.stderr.write("regions lengths do not match\n")
        return (subMatrix, subRegions, regions_no_score)

    @staticmethod
    def coverageFromArray(valuesArray, zones, binSize, avgType):
        try:
            valuesArray[0]
        except IndexError, TypeError:
            sys.stderr.write("values array {}, zones {}\n".format(valuesArray,
                                                                  zones))

        cvgList = []
        start = zones[0][0]
        for zone_start, zone_end, num_bins in zones:
            # the linspace is to get equally spaced positions along the range
            # If the gene is short the sampling regions could overlap,
            # if it is long, the sampling regions would be spaced
            countsList = []

            # this case happens when the downstream or upstream
            # region is set to 0
            if zone_start == zone_end:
                continue

            (posArray, stepSize) = np.linspace(zone_start, zone_end, num_bins,
                                               endpoint=False,
                                               retstep=True)
            stepSize = np.ceil(stepSize)

            for pos in np.floor(posArray):
                indexStart = int(pos - start)
                #indexEnd   = int(indexStart + binSize)
                indexEnd   = int(indexStart + stepSize + 1)
                try:
                    countsList.append(
                        heatmapper.myAverage(valuesArray[indexStart:indexEnd],
                                             avgType))
                except Exception as detail:
                    sys.stderr.write("Exception found. "
                                     "Message: {}\n".format(detail))
            cvgList.append(np.array(countsList))
        return np.concatenate(cvgList)

    @staticmethod
    def changeChromNames(chrom):
        """
        Changes UCSC chromosome names to ensembl chromosome names
        and vice versa.
        TODO: mapping from chromosome names ... e.g. mt, unknown_ ...
        """
        if chrom.startswith('chr'):
            return chrom[3:]
        else:
            return 'chr%s' % chrom

    @staticmethod
    def coverageFromBam(bamfile, chrom, zones, binSize, avgType):
        """
        currently this method is deactivated because is too slow.
        It is preferred to create a coverage bigiwig file from the
        bam file and then run heatmapper.
        """
        if chrom not in bamfile.references:
            chrom = heatmapper.changeChromNames(chrom)
            if chrom not in bamfile.references:
                sys.stderr.write(
                    "Skipping region located at unknown chromosome: {} "
                    "Known chromosomes are: {}\n".format(chrom,
                                                         bamfile.references))
                return None
            else:
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
            valuesArray = np.zeros(end - start)
            for read in bamfile.fetch(chrom, min(0, start), end):
                indexStart = max(read.pos - start, 0)
                indexEnd = min(read.pos - start + read.qlen, end - start)
                valuesArray[indexStart:indexEnd] += 1
        except ValueError:
            sys.stderr.write(
                "Value out of range for region {}s {} {}"
                "\n".format(chrom, start, end))
            return np.array([0])  # return something inocuous

        return heatmapper.coverageFromArray(valuesArray, zones,
                                            binSize, avgType)

    @staticmethod
    def coverageFromBigWig(bigwig, chrom, zones, binSize, avgType,
                           nansAsZeros=False):

        """
        uses bigwig file reader from bx-python
        to query a region define by chrom and zones.
        The output is an array that contains the bigwig
        value per base pair. The summary over bins is
        done in a later step when coverageFromArray is called.
        This method is more reliable than quering the bins
        directly from the bigwig, which should be more efficient.

        By default, any region, even if no chromosome match is found
        on the bigwig file, produces a result. In other words
        no regions are skipped.

        This is useful if several matrices wants to be merged
        or if the sorted BED output of one computeMatrix operation
        needs to be used for other cases
        """

        # intialize values array. The length of the array
        # is the length of the region which is defined
        # by the start of the first zone zones[0][0]
        # to the end of the last zone zones[-1][1]
        valuesArray = np.zeros(zones[-1][1] - zones[0][0])
        if not nansAsZeros:
            valuesArray[:] = np.nan
        try:
            bw_array = bigwig.get_as_array(chrom,
                                           max(0, zones[0][0]),
                                           zones[-1][1])
        except Exception as detail:
                sys.stderr.write("Exception found. Message: "
                                 "{}\n".format(detail))

        if bw_array is None:
            # When bigwig.get_as_array queries a
            # chromosome that is not known
            # it returns None. Ideally, the bigwig should
            # be able to inform the known chromosome names
            # as is the case for bam files, but the
            # bx-python function does not allow access to
            # this info.
            altered_chrom = heatmapper.changeChromNames(chrom)
            bw_array = bigwig.get_as_array(altered_chrom,
                                           max(0, zones[0][0]),
                                           zones[-1][1])
            # test again if with the altered chromosome name
            # the bigwig returns something.
            if bw_array is None:
                sys.stderr.write("Warning: Your chromosome names do "
                                 "not match.\nPlease check that the "
                                 "chromosome names in your BED "
                                 "file\ncorrespond to the names in your "
                                 "bigWig file.\nAn empty line will be "
                                 "added you your heatmap.\nThe offending "
                                 "chromosome name is "
                                 "{}\n\n".format(chrom))

        if bw_array is not None:
            if zones[0][0] < 0:
                valuesArray = np.zeros(zones[-1][1] - zones[0][0])
                valuesArray[:] = np.nan
                valuesArray[abs(zones[0][0]):] = bw_array
            else:
                valuesArray = bw_array

        # replaces nans for zeros
        if nansAsZeros:
            valuesArray[np.isnan(valuesArray)] = 0
        return heatmapper.coverageFromArray(valuesArray, zones,
                                            binSize, avgType)

    @staticmethod
    def myAverage(valuesArray, avgType='mean'):
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

    def matrixFromDict(self, matrixDict, regionsDict, parameters):
        self.regionsDict = regionsDict
        self.matrixDict = matrixDict
        self.parameters = parameters
        self.lengthDict = OrderedDict()
        self.matrixAvgsDict = OrderedDict()

    def readMatrixFile(self, matrix_file, verbose=None,
                       default_group_name='label_1'):
        # reads a bed file containing the position
        # of genomic intervals
        # In case a hash sign '#' is found in the
        # file, this is considered as a delimiter
        # to split the heatmap into groups

        regions = []
        matrix_rows = []
        regionsDict = OrderedDict()
        matrixDict = OrderedDict()
        matrixAvgsDict = OrderedDict()
        regionGroups = [(0, '')]
        parameters = dict()
        totalIntervals = 0
        includedIntervals = 0

        fh = gzip.open(matrix_file)
        for line in fh:
            line = line.strip()
            totalIntervals += 1
            # read the header file containing the parameters
            # used
            if line.startswith("@"):
                # the parameters used are saved using the
                # a key1:value1\tkey2:value2 format
                for p in line[1:].strip().split('\t'):
                    key, value = p.split(":")
                    try:
                        value = int(value)
                    except ValueError:
                        pass
                    if value == 'True':
                        value = True
                    elif value == 'False':
                        value = False
                    elif value == 'None':
                        value = None
                    parameters[key] = value
                continue
            if line.startswith('#'):
                if includedIntervals > 0 and  \
                        includedIntervals - regionGroups[-1][0] > 0:
                    label = line[1:]
                    regionsDict[label] = np.array(regions[:])
                    matrixDict[label] = \
                        np.ma.masked_invalid(np.vstack(matrix_rows))
                    # be default compute the mean average
                    # this will be rewritten in case the user
                    # has choosen to sort the matrix
                    matrixAvgsDict[label] = np.mean(matrixDict[label], axis=1)
                    regions = []
                    matrix_rows = []
                continue
            region = line.split('\t')
            chrom, start, end, name, mean, strand = region[0:6]
            matrix_rows.append(np.fromiter(region[6:], np.float))
            regions.append({'chrom': chrom, 'start': int(start),
                            'end': int(end), 'name': name, 'mean': float(mean),
                            'strand': strand})
            includedIntervals += 1

        if len(regions):
            regionsDict[default_group_name] = np.array(regions)
            matrixDict[default_group_name] = \
                np.ma.masked_invalid(np.vstack(matrix_rows))
            matrixAvgsDict[label] = np.mean(matrixDict[default_group_name],
                                            axis=1)

        self.regionsDict = regionsDict
        self.matrixDict = matrixDict
        self.parameters = parameters
        self.lengthDict = OrderedDict()
        self.matrixAvgsDict = OrderedDict()

        return

    def sortMatrix(self, sort_using='mean', sort_method='no'):
        # sort the matrix using the average of values per row
        self.matrixAvgsDict = OrderedDict()
        self.lengthDict = OrderedDict()
        for label in self.matrixDict.keys():
            if sort_method != 'no':
                if sort_using == 'region_length':
                    matrixAvgs = np.array([x['end'] - x['start']
                                           for x in self.regionsDict[label]])
                    b = self.parameters['upstream'] / \
                        self.parameters['bin size']
                    # for plotting I add the upstream
                    # distance
                    self.lengthDict[label] = \
                        b + (matrixAvgs / self.parameters['bin size'])
                else:
                    matrixAvgs = np.__getattribute__(sort_using)(
                        self.matrixDict[label], axis=1)
                SS = matrixAvgs.argsort()

                if sort_method == 'descend':
                    SS = SS[::-1]
                self.matrixDict[label] = self.matrixDict[label][SS, :]
                self.regionsDict[label] = self.regionsDict[label][SS]
                self.matrixAvgsDict[label] = matrixAvgs[SS]
            try:
                self.lengthDict[label] = self.lengthDict[label][SS]
            except:
                self.lengthDict[label] = None

    def saveMatrix(self, file_name):
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
        fh = gzip.open(file_name, 'wb')
        params_str = '\t'.join(["{}:{}".format(key, value)
                                for key, value in self.parameters.iteritems()])
        fh.write("@" + params_str + "\n")
        for label, regions in self.regionsDict.iteritems():
            j = 0
            for region in regions:
                # this method to join np_array values
                # keeps nans while converting them to strings
                score = 0
                # score = region['mean']
                if self.matrixAvgsDict is not None:
                    if np.ma.is_masked(self.matrixAvgsDict[label][j]):
                        score = 'nan'
                    else:
                        score = np.float(self.matrixAvgsDict[label][j])
                matrix_values = "\t".join(
                    np.char.mod('%f', self.matrixDict[label][j]))
                fh.write(
                    '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        region['chrom'],
                        region['start'],
                        region['end'],
                        region['name'],
                        score,
                        region['strand'],
                        matrix_values))
                j += 1
            fh.write('#{}\n'.format(label))
        fh.close()

    def saveTabulatedValues(self, file_handle):
        bin = range(self.parameters['upstream'] * -1,
                    self.parameters['body'] + self.parameters['downstream'],
                    self.parameters['bin size'])

        avgDict = OrderedDict()
        stdDict = OrderedDict()

        for label, heatmapMatrix in self.matrixDict.iteritems():
            avgDict[label] = heatmapper.matrixAvg(heatmapMatrix, 'mean')
            stdDict[label] = heatmapper.matrixAvg(heatmapMatrix, 'std')

        file_handle.write(
            '#bin No.\t{}\n'.format(" mean\t std\t".join(avgDict.keys())))

        for j in range(0, len(avgDict[avgDict.keys()[0]])):
            file_handle.write('{}\t'.format(bin[j]))
            for label in self.matrixDict.keys():
                file_handle.write(
                    '{}\t{}\t'.format(avgDict[label][j], stdDict[label][j]))
            file_handle.write('\n')

        file_handle.close()

    def reLabelGroups(self, group_labels):
        """
        changes the labels of the groups
        """
        if len(group_labels) != len(self.matrixDict):
            sys.stderr.write(
                "The number of groups does not match the number of labels\n"
                "given. Please define {} different labels.\nThe labels\n"
                "given were: {}\n.".format(len(self.matrixDict), group_labels))
            exit(1)
        elif len(set(group_labels)) != len(group_labels):
            sys.stderr.write(
                "The group labels given contain repeated names. Please\n"
                "give a unique name to each value. The values given are\n"
                "{} \n".format(group_labels))
            exit(1)

        matrixDict = OrderedDict()
        regionsDict = OrderedDict()
        lengthDict = OrderedDict()
        matrixAvgsDict = OrderedDict()
        old_labels = self.matrixDict.keys()
        for index in range(len(group_labels)):
            matrixDict[group_labels[index]] = \
                self.matrixDict[old_labels[index]]
            regionsDict[group_labels[index]] = \
                self.regionsDict[old_labels[index]]
            try:
                lengthDict[group_labels[index]] = \
                    self.lengthDict[old_labels[index]]
            except KeyError:
                pass
            try:
                matrixAvgsDict[group_labels[index]] = \
                    self.matrixAvgsDict[old_labels[index]]
            except KeyError:
                pass

        self.matrixDict = matrixDict
        self.regionsDict = regionsDict
        self.lengthDict = lengthDict
        self.matrixAvgsDict = matrixAvgsDict

    def saveMatrixValues(self, file_name):
        # print a header telling the group names and their length
        fh = open(file_name, 'w')
        info = []
        for label, regions in self.regionsDict.iteritems():
            info.append("{}:{}".format(label, len(regions)))
        fh.write("#{}\n".format("\t".join(info)))
        # add to header the x axis values
        fh.write("#downstream:{}\tupstream:{}\tbody:{}\tbin size:{}\n".format(
                 self.parameters['downstream'],
                 self.parameters['upstream'],
                 self.parameters['body'],
                 self.parameters['bin size']))

        fh.close()
        # reopen again using append mode
        fh = open(file_name, 'a')
        for key in self.matrixDict:
            np.savetxt(fh, self.matrixDict[key], fmt="%.3g")
        fh.close()

    def saveBED(self, file_handle):
        for label, regions in self.regionsDict.iteritems():
            cluster = label.find('cluster') != -1
            j = 0
            for region in regions:
                score = 0
                if self.matrixAvgsDict is not None:
                    try:
                        if self.matrixAvgsDict[label][j] is np.ma.masked:
                            score = 'nan'
                        else:
                            score = np.float(self.matrixAvgsDict[label][j])
                    except KeyError:
                        pass
                # If the label is from clustering, we add an additional column
                # TODO: we need to decide if you want to go with an additional column
                # also for other labels. Comment lines are not in the BED specification
                # and can cause additional errors.
                if cluster:
                    file_handle.write(
                        '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            region['chrom'],
                            region['start'],
                            region['end'],
                            region['name'],
                            score,
                            region['strand'],
                            label))
                else:
                    file_handle.write(
                        '{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            region['chrom'],
                            region['start'],
                            region['end'],
                            region['name'],
                            score,
                            region['strand']))
                j += 1
            if not cluster:
                file_handle.write('#{}\n'.format(label))
        file_handle.close()

    @staticmethod
    def matrixAvg(matrix, avgType='mean'):
        matrix = np.ma.masked_invalid(matrix)
        return np.__getattribute__(avgType)(matrix, axis=0)

    @staticmethod
    def filterGenomicIntervalFile(file_handle):
        """
        Filter track lines out of a GenomicIntervalFile, normally from UCSC.
        Return an iterator over the lines of file_handle.
        """
        for line in file_handle:
            if line.startswith('browser') or line.startswith('track'):
                continue
            yield line

    @staticmethod
    def getRegionsAndGroups(regions_file, onlyMultiplesOf=1,
                            default_group_name='genes',
                            verbose=None):
        # reads a bed file containing the position
        # of genomic intervals
        # In case is hash sign '#' is found in the
        # file, this is considered as a delimiter
        # to split the heatmap into groups

        regions = []
        regionsDict = OrderedDict()
        regionGroups = [(0, '')]

        prevInterval = None
        duplicates = 0
        totalIntervals = 0
        includedIntervals = 0
        for ginterval in GenomicIntervalReader(
            heatmapper.filterGenomicIntervalFile(regions_file),
            fix_strand=True):

            totalIntervals += 1
            if ginterval.__str__().startswith('#'):
                if includedIntervals > 1 and  \
                        includedIntervals - regionGroups[-1][0] > 1:
                    label = ginterval.__str__()[1:]
                    newLabel = label
                    if label in regionsDict.keys():
                       # loop to find a unique label name
                        i = 0
                        while True:
                            i += 1
                            newLabel = label + "_r" + str(i)
                            if newLabel not in regionsDict.keys():
                                break

                    regionsDict[newLabel] = np.array(regions[:])
                    regions = []
                continue
            # if the list of regions is to big, only
            # consider a fraction of the data
            if totalIntervals % onlyMultiplesOf != 0:
                continue
            # check for regions that have the same position as the previous.
            # This assumes that the regions file given is sorted
            if prevInterval and prevInterval.chrom == ginterval.chrom and \
                    prevInterval.start == ginterval.start and \
                    prevInterval.end == ginterval.end and \
                    prevInterval.strand == ginterval.strand:
                if verbose:
                    try:
                        genename = ginterval.fields[3]
                    except:
                        genename = ''
                    sys.stderr.write("*Warning* Duplicated region: "
                                     "{} {}:{}-{}.\n".format(
                            genename,
                            ginterval.chrom, ginterval.start,
                            ginterval.end))
                duplicates += 1
            else:
                prevInterval = ginterval

            regions.append(heatmapper.ginterval2dict(ginterval))
            includedIntervals += 1

        # in case we reach the end of the file
        # without encountering a hash,
        # a default name is given to regions
        if len(regions):
            regionsDict[default_group_name] = np.array(regions)

        if verbose and duplicates > 0:
            sys.stderr.write(
                "{} ({:.2f}) regions covering the exact same interval\n"
                "were found".format(duplicates,
                                    float(duplicates) * 100 / totalIntervals))

        return regionsDict

    @staticmethod
    def ginterval2dict(genomicInterval):
        """
        transforms a genomic interval from bx python
        into a dictionary
        """
        region = {'chrom': genomicInterval.chrom,
                  'start': genomicInterval.start,
                  'end': genomicInterval.end,
                  'strand': genomicInterval.strand}
        try:
            region['name'] = genomicInterval.fields[3]
        except IndexError:
            region['name'] = "No name"
        return region

    @staticmethod
    def hmcluster(matrix, k, method='kmeans'):

        matrix = np.asarray(matrix)
        # replace nans for 0 otherwise kmeans produces a weird behaviour
        matrix[np.isnan(matrix)] = 0

        if method == 'kmeans':
            from scipy.cluster.vq import vq, kmeans

            centroids, _ = kmeans(matrix, k)
            # order the centroids in an attempt to
            # get the same cluster order
            order = np.argsort(centroids.mean(axis=1))
            idx,_ = vq(matrix, centroids[order,])

        if method == 'hierarchical':
            from scipy.cluster.hierarchy import fcluster, linkage
            Z = linkage(matrix, method='ward')
            idx = fcluster(Z, k, criterion='maxclust')

        return idx
