import sys


class ReadBed(object):

    """
    Reads a bed file. Based on the number of fields
    it tries to guess the type of bed file used. Current options
    are bed3, bed6 and bed12

    Examples
    --------
    bed = readBed(open("file.bed", 'r'))
    for interval in bed:
    ... print bed['start']

    """

    def __init__(self, file_handle):
        """
        Parameters
        ----------
        file_name : str
            name of the bed file
        """

        self.file_type = None
        self.fields = ['chrom', 'start', 'end',
                       'name', 'score', 'strand',
                       'thickStart', 'thickEnd',
                       'itemRgb', 'blockCount',
                       'blockSizes', 'blockStarts']

        self.file_handle = file_handle

    def __iter__(self):
        return self

    def get_no_comment_line(self):
        """
        Skips comment lines starting with '#'
        "track" or "browser" in the bed files
        """
        line = self.file_handle.readline()
        if line.startswith("#") or line.startswith("track") or \
           line.startswith("browser"):
            line = self.get_no_comment_line()

        return line

    def guess_file_type(self, line_values):
        """try to guess type of bed file by counting the fields
        """
        if len(line_values) == 3:
            self.file_type = 'bed3'
        elif len(line_values) == 6:
            self.file_type = 'bed6'
        elif len(line_values) == 12:
            self.file_type = 'bed12'
        elif len(line_values) > 6:
            # assume bed6
            self.file_type = 'bed6'
            sys.stderr.write("Number of fields in BED file is not standard. Assuming bed6\n")
        else:
            # assume bed3
            self.file_type = 'bed3'
            sys.stderr.write("Number of fields in BED file is not standard. Assuming bed3\n")
        return self.file_type

    def next(self):
        """
        :return: bedInterval object
        """
        line = self.file_handle.next()

        # skip empty lines
        while True:
            if line.strip() == '':
                line = self.file_handle.next()
            else:
                break
        if line.startswith("#") or line.startswith("track") or \
           line.startswith("browser"):
            # heatmapper identifies clusters by looking at
            # the '#' values
            return BedInterval(dict(), line.strip())

        return BedInterval(self.get_values_from_line(line), line)

    def get_values_from_line(self, bed_line):
        """
        Processes each bed line from a bed file
        and casts the values

        """

        line_data = bed_line.strip().split("\t")

        line_values = []
        for idx, r in enumerate(line_data):
            # first field is always chromosome name
            # or scaffold name and should be cast as a string
            # same for field 3 (name) and field 6 (strand)
            if idx == 0 or idx == 3 or idx == 5:
                line_values.append(r)
            elif idx == 1 or idx == 2:
                # start and end fields must be integers
                try:
                    line_values.append(int(r))
                except ValueError:
                    sys.stderr.write("Value: {} in field {} is not an integer".format(r, idx + 1))
                    return dict()
            else:
                tmp = r
                try:
                    tmp = float(r)
                except (ValueError, TypeError):
                    tmp = r
                line_values.append(tmp)
        if self.file_type is None:
            self.file_type = self.guess_file_type(line_values)

        if self.file_type == 'bed3':
            line_values = line_values[0:3]
            # in case of a bed3, the id, score and strand
            # values are added as ".", 0, "." respectively
            line_values.extend([".", 0, "."])
        elif self.file_type == 'bed6':
                line_values = line_values[0:6]

        fields_dict = {}
        try:
            fields_dict = dict([(self.fields[index], value) for index, value in enumerate(line_values)])
        except:
            import ipdb
            ipdb.set_trace()

        return fields_dict


class BedInterval(object):

    """
    simple place holder for the line data of a bed file
    """

    def __init__(self, field_data_dict, line):
        """
        Parameters
        ----------
        field_data_dict : dict
            A dictionary of field values. E.g. {'chrom': 'chr1', 'start': 0 } etc.
        line : str
            string of the file line.

        >>> bed = BedInterval({'chrom': '1', 'start': 1, 'end': 10}, "")
        >>> bed.start
        1
        >>> bed['start']
        1
        """
        self.fields = {}
        for key, value in field_data_dict.iteritems():
            setattr(self, key, value)
            self.fields[key] = value

        self.line = line

    def __getitem__(self, item):
        return self.fields[item]
