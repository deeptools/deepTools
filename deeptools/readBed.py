import sys


class ReadBed(object):
    """
    Reads a bed file. Based on the number of fields
    it tries to guess the type of bed file used. Current options
    are bed3, bed6 and bed12

    Example:
    bed = readBed(open("file.bed", 'r'))
    for interval in bed:
        print bed['start']

    """

    def __init__(self, file_handle):

        """
        :param file_name: name of the bed file
        :return:
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
        :return:
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
            sys.stderr.write("Number of fields in BED file is not standard. Assuming bed6")
        else:
            # assume bed3
            self.file_type = 'bed3'
            sys.stderr.write("Number of fields in BED file is not standard. Assuming bed3")

    def next(self):
        """
        :return: bedInterval object
        """
        line = self.file_handle.next()

        # skip empty lines
        while 1:
            if line.strip() == '':
                import ipdb;ipdb.set_trace()
                line = self.file_handle.next()
            else:
                break
        if line.startswith("#") or line.startswith("track") or \
           line.startswith("browser"):
            return BedInterval(dict(), line.strip())

        line_data = line.strip().split("\t")

        line_values = []
        for r in line_data:
            if r.isdigit():
                line_values.append(int(r))
            else:
                tmp = r
                try:
                    tmp = float(r)
                except ValueError, TypeError:
                    tmp = r
                line_values.append(tmp)
        if self.file_type is None:
            self.guess_file_type(line_values)

        fields_dict = dict([(self.fields[index], value) for index, value in enumerate(line_values)])
        return BedInterval(fields_dict, line)


class BedInterval(object):
    """
    simple place holder for the line data of a bed file
    """

    def __init__(self, field_data_dict, line):
        """

        :param field_data_dict: A dictionary of field values. E.g. {'chrom': 'chr1', 'start': 0 } etc.
        :param line: string of the file line.
        """
        self.fields = {}
        for key, value in field_data_dict.iteritems():
            setattr(self, key, value)
            self.fields[key] = value

        self.line  = line

    def __getitem__(self, item):
        return self.fields[item]
