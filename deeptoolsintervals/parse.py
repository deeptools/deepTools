#!/usr/bin/env python

from deeptoolsintervals import tree
import sys
import gzip
try:
    import bz2
    supportsBZ2 = True
except:
    supportsBZ2 = False
import os.path
import csv


def getNext(fp):
    """
    Sometimes we need to decode, sometimes not
    """
    line = fp.readline()
    if isinstance(line, str):
        return line
    return line.decode('ascii')


def seemsLikeGTF(cols):
    """
    Does a line look like it could be from a GTF file? Column contents must be:

    3: int
    4: int
    5: '.' or float
    6: '+', '-', or '.'
    7: 0, 1 or 2
    8: matches the attribute regular expression
    """
    try:
        int(cols[3])
        int(cols[4])
        if cols[5] != '.':
            float(cols[5])
        cols[6] in ['+', '-', '.']
        if cols[7] != '.':
            int(cols[7]) in [0, 1, 2]
        s = next(csv.reader([cols[8]], delimiter=' '))
        assert("gene_id" in s)
        assert(s[-1] != "gene_id")
        return True
    except:
        return False


def findRandomLabel(labels, name):
    """
    Because some people are too clever by half, ensure that group labels are unique...
    """
    if name not in labels:
        return name

    # This is what the heatmapper.py did to ensure unique names
    i = 0
    while True:
        i += 1
        nameTry = name + "_r" + str(i)
        if nameTry not in labels:
            return nameTry


def parseExonBounds(start, end, n, sizes, offsets):
    """
    Parse the last 2 columns of a BED12 file and return a list of tuples with
    (exon start, exon end) entries.

    If the line is malformed, issue a warning and return (start, end)
    """
    offsets = offsets.strip(",").split(",")
    sizes = sizes.strip(",").split(",")
    offsets = offsets[0:n]
    sizes = sizes[0:n]
    try:
        starts = [start + int(x) for x in offsets]
        ends = [start + int(x) + int(y) for x, y in zip(offsets, sizes)]
    except:
        sys.stderr.write("Warning: Received an invalid exon offset ({0}) or size ({1}), using the entry bounds instead ({2}-{3})\n".format(offsets, sizes, start, end))
        return [(start, end)]

    if len(offsets) < n or len(sizes) < n:
        sys.stderr.write("Warning: There were too few exon start/end offsets ({0}) or sizes ({1}), using the entry bounds instead ({2}-{3})\n".format(offsets, sizes, start, end))
        return [(start, end)]

    return [(x, y) for x, y in zip(starts, ends)]


def openPossiblyCompressed(fname):
    """
    A wrapper to open gzip/bzip/uncompressed files
    """
    with open(fname, "rb") as f:
        first3 = bytes(f.read(3))
    if first3 == b"\x1f\x8b\x08":
        return gzip.open(fname, "rb")
    elif first3 == b"\x42\x5a\x68" and supportsBZ2:
        return bz2.BZ2File(fname, "rb")
    else:
        return open(fname)


def getLabel(line):
    """
    Split by tabs and return the index of "deepTools_group" (or None)
    """
    cols = line.strip().split("\t")
    if "deepTools_group" in cols:
        return cols.index("deepTools_group")
    return None


class GTF(object):
    """
    A class to hold an interval tree and its associated functions

    >>> from deeptoolsintervals import parse
    >>> from os.path import dirname
    >>> gtf = parse.GTF("{0}/test/GRCh38.84.gtf.gz".format(dirname(parse.__file__)), keepExons=True)
    >>> gtf.findOverlaps("1", 1, 20000)
    [(11868, 14409, 'ENST00000456328', 'group 1', [(11868, 12227), (12612, 12721), (13220, 14409)], '.'), (12009, 13670, 'ENST00000450305', 'group 1', [(12009, 12057), (12178, 12227), (12612, 12697), (12974, 13052), (13220, 13374), (13452, 13670)], '.'), (14403, 29570, 'ENST00000488147', 'group 1', [(14403, 14501), (15004, 15038), (15795, 15947), (16606, 16765), (16857, 17055), (17232, 17368), (17605, 17742), (17914, 18061), (18267, 18366), (24737, 24891), (29533, 29570)], '.'), (17368, 17436, 'ENST00000619216', 'group 2', [(17368, 17436)], '.')]
    >>> gtf = parse.GTF("{0}/test/GRCh38.84.gtf.gz".format(dirname(parse.__file__)))
    >>> gtf.findOverlaps("1", 1, 20000)
    [(11868, 14409, 'ENST00000456328', 'group 1', [(11868, 14409)], '.'), (12009, 13670, 'ENST00000450305', 'group 1', [(12009, 13670)], '.'), (14403, 29570, 'ENST00000488147', 'group 1', [(14403, 29570)], '.'), (17368, 17436, 'ENST00000619216', 'group 2', [(17368, 17436)], '.')]
    >>> gtf.findOverlaps("1", 12000, 20000, trimOverlap=True)
    [(12009, 13670, 'ENST00000450305', 'group 1', [(12009, 13670)], '.'), (14403, 29570, 'ENST00000488147', 'group 1', [(14403, 29570)], '.'), (17368, 17436, 'ENST00000619216', 'group 2', [(17368, 17436)], '.')]
    >>> gtf.findOverlaps("1", 1, 20000, numericGroups=True, includeStrand=True)
    [(11868, 14409, 'ENST00000456328', 0, [(11868, 14409)], '+', '.'), (12009, 13670, 'ENST00000450305', 0, [(12009, 13670)], '+', '.'), (14403, 29570, 'ENST00000488147', 0, [(14403, 29570)], '-', '.'), (17368, 17436, 'ENST00000619216', 1, [(17368, 17436)], '-', '.')]
    """

    def firstNonComment(self, fp):
        """
        Skip lines at the beginning of a file starting with #, browser, or track.
        Returns a tuple of the first non-comment line and the column holding the group label (if it exists)
        """
        line = getNext(fp)
        labelColumn = None
        try:
            while line.startswith("#") or line.startswith('track') or line.startswith('browser'):
                if labelColumn is None:
                    labelColumn = getLabel(line)
                line = getNext(fp)
        except:
            sys.stderr.write("Warning, {0} was empty\n".format(self.filename))
            return None
        return line, labelColumn

    def inferType(self, fp, line, labelColumn=None):
        """
        Attempt to infer a file type from a single line. This is largely based on the number of columns plus looking for "gene_id".
        """
        subtract = 0
        if labelColumn is not None:
            subtract = 1
        cols = line.split("\t")
        if len(cols) - subtract < 3:
            raise RuntimeError('{0} does not seem to be a recognized file type!'.format(self.filename))
        elif len(cols) - subtract == 3:
            return 'BED3'
        elif len(cols) - subtract < 6:
            if self.verbose:
                sys.stderr.write("Warning, {0} has an abnormal number of fields. Assuming BED3 format.\n".format(self.filename))
            return 'BED3'
        elif len(cols) - subtract == 6:
            return 'BED6'
        elif len(cols) and seemsLikeGTF(cols):
            return 'GTF'
        elif len(cols) - subtract == 12:
            return 'BED12'
        elif len(cols) - subtract < 12:
            if self.verbose:
                sys.stderr.write("Warning, {0} has an abnormal format. Assuming BED6 format.\n".format(self.filename))
            return 'BED6'
        else:
            if self.verbose:
                sys.stderr.write("Warning, {0} has an abnormal format. Assuming BED12 format.\n".format(self.filename))
            return 'BED12'

    def mungeChromosome(self, chrom, append=True):
        """
        Return the chromosome name, possibly munged to match one already found in the chromosome dictionary
        """
        if chrom in self.chroms:
            return chrom

        # chrM <-> MT and chr1 <-> 1 conversions
        if chrom == "MT" and "chrM" in self.chroms:
            chrom = "chrM"
        elif chrom == "chrM" and "MT" in self.chroms:
            chrom = "MT"
        elif chrom.startswith("chr") and len(chrom) > 3 and chrom[3:] in self.chroms:
            chrom = chrom[3:]
        elif "chr" + chrom in self.chroms:
            chrom = "chr" + chrom

        if append:
            self.chroms.append(chrom)

        return chrom

    def parseBEDcore(self, line, ncols):
        """
        Returns True if the entry was added, otherwise False

        >>> from deeptoolsintervals import parse
        >>> from os.path import dirname
        >>> gtf = parse.GTF("{0}/test/GRCh38.84.bed12.bz2".format(dirname(parse.__file__)), keepExons=True, labels=["foo"])
        >>> gtf.findOverlaps("1", 1, 20000)
        [(11868, 14409, 'ENST00000456328.2', 'foo', [(11868, 12227), (12612, 12721), (13220, 14409)], 0.0), (12009, 13670, 'ENST00000450305.2', 'foo', [(12009, 12057), (12178, 12227), (12612, 12697), (12974, 13052), (13220, 13374), (13452, 13670)], 0.0), (14403, 29570, 'ENST00000488147.1', 'foo', [(14403, 14501), (15004, 15038), (15795, 15947), (16606, 16765), (16857, 17055), (17232, 17368), (17605, 17742), (17914, 18061), (18267, 18366), (24737, 24891), (29533, 29570)], 0.0), (17368, 17436, 'ENST00000619216.1', 'foo', [(17368, 17436)], 0.0)]
        """

        strand = 3
        cols = line.split("\t")
        name = "{0}:{1}-{2}".format(cols[0], cols[1], cols[2])

        if int(cols[1]) < 0:
            cols[1] = 0

        if int(cols[1]) >= int(cols[2]):
            sys.stderr.write("Warning: {0}:{1}-{2} is an invalid BED interval! Ignoring it.\n".format(cols[0], cols[1], cols[2]))
            return

        # BED6/BED12: set name and strand
        score = '.'
        if ncols > 3:
            name = cols[3]
            if cols[5] == '+':
                strand = 0
            elif cols[5] == '-':
                strand = 1
            score = cols[4]

        # Ensure that the name is unique
        name = findRandomLabel(self.exons, name)
        self.tree.addEntry(self.mungeChromosome(cols[0]), int(cols[1]), int(cols[2]), name, strand, self.labelIdx, score)
        if ncols != 12 or self.keepExons is False:
            self.exons[name] = [(int(cols[1]), int(cols[2]))]
        else:
            assert(len(cols) == 12)
            self.exons[name] = parseExonBounds(int(cols[1]), int(cols[2]), int(cols[9]), cols[10], cols[11])

    def parseBED(self, fp, line, ncols=3, labelColumn=None):
        """
        parse a BED file. The default group label is the file name.

        fp:    A python file pointer
        line:  The first line
        ncols: The number of columns to care about


        >>> from deeptoolsintervals import parse
        >>> from os.path import dirname, basename
        >>> gtf = parse.GTF("{0}/test/GRCh38.84.bed6".format(dirname(parse.__file__)), keepExons=True)
        >>> gtf.findOverlaps("1", 1, 20000)
        [(11868, 14409, 'ENST00000456328.2', 'group 1', [(11868, 14409)], 0.0), (12009, 13670, 'ENST00000450305.2', 'group 1', [(12009, 13670)], 0.0), (14403, 29570, 'ENST00000488147.1', 'group 1', [(14403, 29570)], 0.0), (17368, 17436, 'ENST00000619216.1', 'group 1', [(17368, 17436)], 0.0)]
        >>> gtf = parse.GTF("{0}/test/GRCh38.84.bed".format(dirname(parse.__file__)), keepExons=True, labels=["foo", "bar", "quux", "sniggly"])
        >>> gtf.findOverlaps("1", 1, 20000)
        [(11868, 14409, '1:11868-14409', 'foo', [(11868, 14409)], '.'), (12009, 13670, '1:12009-13670', 'foo', [(12009, 13670)], '.'), (14403, 29570, '1:14403-29570', 'foo', [(14403, 29570)], '.'), (17368, 17436, '1:17368-17436', 'foo', [(17368, 17436)], '.')]

        Test having a header in one file, but not another:

        >>> gtf = parse.GTF(["{0}/test/GRCh38.84.labels.bed".format(dirname(parse.__file__)), "{0}/test/GRCh38.84.bed2".format(dirname(parse.__file__))])
        >>> overlaps = gtf.findOverlaps("1", 1, 30000000)
        >>> labels = dict()
        >>> for o in overlaps:
        ...     if basename(o[3]) not in labels:
        ...         labels[basename(o[3])] = 0
        ...     labels[basename(o[3])] += 1
        >>> assert(labels['group 1'] == 4)
        >>> assert(labels['group 2'] == 9)
        >>> assert(labels['group 3'] == 7)
        >>> assert(labels['group 4'] == 1)
        >>> assert(labels['group 1_r1'] == 4)
        >>> assert(labels['group2'] == 9)
        >>> assert(labels['group 3'] == 7)
        >>> assert(labels['GRCh38.84.bed2'] == 1)
        >>> gtf = parse.GTF(["{0}/test/GRCh38.84.bed2".format(dirname(parse.__file__)), "{0}/test/GRCh38.84.labels.bed".format(dirname(parse.__file__))])
        >>> overlaps = gtf.findOverlaps("1", 1, 30000000)
        >>> labels = dict()
        >>> for o in overlaps:
        ...     if basename(o[3]) not in labels:
        ...         labels[basename(o[3])] = 0
        ...     labels[basename(o[3])] += 1
        >>> assert(labels['group 1'] == 8)
        >>> assert(labels['group 2'] == 9)
        >>> assert(labels['group 3'] == 14)
        >>> assert(labels['group 4'] == 1)
        >>> assert(labels['group2'] == 9)
        >>> assert(labels['GRCh38.84.bed2'] == 1)
        """
        groupLabelsFound = 0
        groupEntries = 0

        # Handle the first line
        if labelColumn is not None:
            cols = line.split("\t")
            label = cols.pop(labelColumn)
            line = "\t".join(cols)
            if label in self.labels:
                self.labelIdx = self.labels.index(label)
            else:
                self.labels.append(label)
                self.labelIdx = len(self.labels) - 1
        self.parseBEDcore(line, ncols)
        groupEntries = 1

        # iterate over the remaining lines
        for line in fp:
            if not isinstance(line, str):
                line = line.decode('ascii')
            line = line.strip()
            if len(line) == 0:
                # Apparently this happens, some people seem to like trying to break things
                continue

            if line.startswith("#") and labelColumn is None:
                # If there was a previous group AND it had no entries then remove it
                if groupLabelsFound > 0:
                    if groupEntries == 0:
                        sys.stderr.write("Warning, the '{0}' group had no valid entries! Removing it.\n".format(self.labels[self.labelIdx]))
                        del self.labels[-1]
                        groupLabelsFound -= 1
                        self.labelIdx -= 1

                label = line[1:].strip()
                if len(label):
                    # Guard against duplicate group labels
                    self.labels.append(findRandomLabel(self.labels, label))
                else:
                    # I'm sure someone will try an empty label...
                    self.labels.append(findRandomLabel(self.labels, os.path.basename(self.filename)))
                self.labelIdx += 1
                groupLabelsFound += 1
                groupEntries = 0
            elif line.startswith("#") and labelColumn is not None:
                continue
            else:
                if labelColumn is not None:
                    cols = line.split("\t")
                    label = cols.pop(labelColumn)
                    line = "\t".join(cols)
                    if label in self.labels:
                        self.labelIdx = self.labels.index(label)
                    else:
                        self.labels.append(label)
                        self.labelIdx = len(self.labels) - 1
                self.parseBEDcore(line, ncols)
                if labelColumn is None:
                    groupEntries += 1

        if groupEntries > 0 and labelColumn is None:
            if self.defaultGroup is not None:
                self.labels.append(findRandomLabel(self.labels, self.defaultGroup))
            else:
                self.labels.append(findRandomLabel(self.labels, os.path.basename(self.filename)))

        # Reset self.labelIdx
        self.labelIdx = len(self.labels)

    def parseGTFtranscript(self, cols, label):
        """
        Parse and add a transcript entry
        """
        if int(cols[3]) - 1 < 0:
            sys.stderr.write("Warning: Invalid start in '{0}', skipping\n".format("\t".join(cols)))
            return

        if len(cols) < 9:
            sys.stderr.write("Warning: non-GTF line encountered! {0}\n".format("\t".join(cols)))
            return

        s = next(csv.reader([cols[8]], delimiter=' '))
        if "deepTools_group" in s and s[-1] != "deepTools_group":
            label = s[s.index("deepTools_group") + 1].rstrip(";")
        elif self.defaultGroup is not None:
            label = self.defaultGroup

        if self.transcript_id_designator not in s or s[-1] == self.transcript_id_designator:
            sys.stderr.write("Warning: {0} is malformed!\n".format("\t".join(cols)))
            return

        name = s[s.index(self.transcript_id_designator) + 1].rstrip(";")
        if name in self.exons:
            sys.stderr.write("Warning: {0} occurs more than once! Only using the first instance.\n".format(name))
            self.transcriptIDduplicated.append(name)
            return

        if int(cols[3]) > int(cols[4]) or int(cols[3]) < 1:
            sys.stderr.write("Warning: {0}:{1}-{2} is an invalid GTF interval! Ignoring it.\n".format(cols[0], cols[3], cols[4]))
            return

        strand = 3
        if cols[6] == '+':
            strand = 0
        elif cols[6] == '-':
            strand = 1

        score = cols[5]

        # Get the label index
        if label not in self.labels:
            self.labels.append(label)
        self.labelIdx = self.labels.index(label)

        chrom = self.mungeChromosome(cols[0])
        self.tree.addEntry(chrom, int(cols[3]) - 1, int(cols[4]), name, strand, self.labelIdx, score)

        # Exon bounds placeholder
        self.exons[name] = []

    def parseGTFexon(self, cols):
        """
        Parse an exon entry and add it to the transcript hash
        """
        if int(cols[3]) - 1 < 0:
            sys.stderr.write("Warning: Invalid start in '{0}', skipping\n".format("\t".join(cols)))
            return

        s = next(csv.reader([cols[8]], delimiter=' '))
        if self.transcript_id_designator not in s or s[-1] == self.transcript_id_designator:
            sys.stderr.write("Warning: {0} is malformed!\n".format("\t".join(cols)))
            return

        name = s[s.index(self.transcript_id_designator) + 1].rstrip(";")
        if name in self.transcriptIDduplicated:
            return
        if name not in self.exons:
            self.exons[name] = []

        self.exons[name].append((int(cols[3]) - 1, int(cols[4])))

    def parseGTF(self, fp, line):
        """
        parse a GTF file. Note that a single label will be used for every entry
        in a file that isn't explicitly labeled with a deepTools_group
        key:values pair in the last column

        fp:   A python file pointer
        line: The first non-comment line

        >>> from deeptoolsintervals import parse
        >>> from os.path import dirname, basename
        >>> gtf = parse.GTF(["{0}/test/GRCh38.84.gtf.gz".format(dirname(parse.__file__)), "{0}/test/GRCh38.84.2.gtf.gz".format(dirname(parse.__file__))], keepExons=True)
        >>> overlaps = gtf.findOverlaps("1", 1, 20000000)
        >>> labels = dict()
        >>> for o in overlaps:
        ...     if basename(o[3]) not in labels:
        ...         labels[basename(o[3])] = 0
        ...     labels[basename(o[3])] += 1
        >>> assert(labels['GRCh38.84.gtf.gz'] == 17)
        >>> assert(labels['GRCh38.84.2.gtf.gz'] == 6)
        >>> assert(labels['group 1'] == 5)
        >>> assert(labels['group 2'] == 3)

        Test GTF and a BED file
        >>> gtf = parse.GTF(["{0}/test/GRCh38.84.gtf.gz".format(dirname(parse.__file__)), "{0}/test/GRCh38.84.bed".format(dirname(parse.__file__))])
        >>> overlaps = gtf.findOverlaps("1", 1, 20000000)
        >>> labels = dict()
        >>> for o in overlaps:
        ...     if basename(o[3]) not in labels:
        ...         labels[basename(o[3])] = 0
        ...     labels[basename(o[3])] += 1
        >>> assert(labels['GRCh38.84.gtf.gz'] == 17)
        >>> assert(labels['group 1'] == 3)
        >>> assert(labels['group 2'] == 1)
        >>> assert(labels['group 1_r1'] == 4)
        >>> assert(labels['group2'] == 9)
        >>> assert(labels['group 3'] == 7)
        >>> assert(labels['GRCh38.84.bed'] == 1)
        >>> gtf = parse.GTF(["{0}/test/GRCh38.84.bed".format(dirname(parse.__file__)), "{0}/test/GRCh38.84.gtf.gz".format(dirname(parse.__file__))])
        >>> overlaps = gtf.findOverlaps("1", 1, 20000000)
        >>> labels = dict()
        >>> for o in overlaps:
        ...     if basename(o[3]) not in labels:
        ...         labels[basename(o[3])] = 0
        ...     labels[basename(o[3])] += 1
        >>> assert(labels['GRCh38.84.gtf.gz'] == 17)
        >>> assert(labels['group 1'] == 7)
        >>> assert(labels['group 2'] == 1)
        >>> assert(labels['group2'] == 9)
        >>> assert(labels['group 3'] == 7)
        >>> assert(labels['GRCh38.84.bed'] == 1)
        """
        file_label = findRandomLabel(self.labels, os.path.basename(self.filename))

        # Handle the first line
        cols = line.split("\t")
        if cols[2].lower() == self.transcriptID.lower():
            self.parseGTFtranscript(cols, file_label)
        elif cols[2].lower() == self.exonID.lower():
            self.parseGTFexon(cols)

        # Handle the remaining lines
        for line in fp:
            if not isinstance(line, str):
                line = line.decode('ascii')
            if not line.startswith('#'):
                cols = line.split("\t")
                if len(cols) == 0:
                    continue

                if cols[2].lower() == self.transcriptID.lower():
                    self.parseGTFtranscript(cols, file_label)
                elif cols[2].lower() == self.exonID.lower() and self.keepExons is True:
                    self.parseGTFexon(cols)

        # Reset self.labelIdx
        self.labelIdx = len(self.labels)

    def __init__(self, fnames, exonID="exon", transcriptID="transcript", keepExons=False, labels=[], transcript_id_designator="transcript_id", defaultGroup=None, verbose=False):
        """
        Driver function to actually parse files. The steps are as follows:

        1) skip to the first non-comment line
        2) Infer the type from that
        3) Call a type-specific processing function accordingly
          * These call the underlying C code for storage
          * These handle chromsome name conversions (python-level)
          * These handle labels (python-level, with a C-level numeric attribute)
        4) Sanity checking (do the number of labels make sense?)

        Required inputs are as follows:

        fnames:	A list of (possibly compressed with gzip or bzip2) GTF or BED files.

        Optional input is:

        exonID:	      For GTF files, the feature column (column 3) label for
                      exons, or whatever else should be stored as exons. The
                      default is 'exon', though one could use 'CDS' instead.
        transcriptID: As above, but for transcripts. The default is
                      'transcript_id'.
        keepExons:    For BED12 and GTF files, exons are ignored by default.
        labels:       A list of group labels.
        transcript_id_designator: For gtf files, this is the key used in a searching for the transcript ID.
                      If one sets transcriptID to 'gene', then
                      transcript_id_designator would need to be changed to
                      'gene_id' or 'gene_name' to extract the gene ID/name from
                      the attributes.
        defaultGroup: The default group name. If None, the file name is used.
        verbose:      Whether to produce warning messages (default: False)
        """
        self.fname = []
        self.filename = ""
        self.chroms = []
        self.exons = {}
        self.labels = []
        self.transcriptIDduplicated = []
        self.tree = tree.initTree()
        self.labelIdx = 0
        self.transcript_id_designator = transcript_id_designator
        self.exonID = exonID
        self.transcriptID = transcriptID
        self.keepExons = keepExons
        self.defaultGroup = defaultGroup
        self.verbose = verbose

        if labels != []:
            self.already_input_labels = True

        if not isinstance(fnames, list):
            fnames = [fnames]

        # Load the files
        for fname in fnames:
            self.filename = fname
            fp = openPossiblyCompressed(fname)
            line, labelColumn = self.firstNonComment(fp)
            if line is None:
                # This will only ever happen if a file is empty or just has a header/comment
                continue
            line = line.strip()

            ftype = self.inferType(fp, line, labelColumn)
            if ftype == 'GTF':
                self.parseGTF(fp, line)
            elif ftype == 'BED3':
                self.parseBED(fp, line, 3, labelColumn)
            elif ftype == 'BED6':
                self.parseBED(fp, line, 6, labelColumn)
            else:
                self.parseBED(fp, line, 12, labelColumn)
            fp.close()

        # Sanity check
        if self.tree.countEntries() == 0:
            raise RuntimeError("None of the input BED/GTF files had valid regions")

        # Replace labels
        if len(labels) > 0:
            if len(labels) != len(self.labels):
                raise RuntimeError("The number of labels found ({0}) does not match the number input ({1})!".format(self.labels, labels))
            else:
                self.labels = labels

        # vine -> tree
        self.tree.finish()

    # findOverlaps()
    def findOverlaps(self, chrom, start, end, strand=".", matchType=0, strandType=0, trimOverlap=False, numericGroups=False, includeStrand=False):
        """
        Given a chromosome and start/end coordinates with an optional strand,
        return a list of tuples comprised of:

         * start
         * end
         * name
         * label
         * [(exon start, exon end), ...]
         * strand (optional)

        If there are no overlaps, return None. This function allows stranded
        searching, though the default is to ignore strand!

        The non-obvious options are defined in gtf.h:
          matchType:  0, GTF_MATCH_ANY
                      1, GTF_MATCH_EXACT
                      2, GTF_MATCH_CONTAIN
                      3, GTF_MATCH_WITHIN
                      4, GTF_MATCH_START
                      5, GTF_MATCH_END

          strandType: 0, GTF_IGNORE_STRAND
                      1, GTF_SAME_STRAND
                      2, GTF_OPPOSITE_STRAND
                      3, GTF_EXACT_SAME_STRAND

        trimOverlap:   If true, this removes overlaps from the 5' end that
                       extend beyond the range requested. This is useful in
                       cases where a function calling this does is first
                       dividing the genome into large bins. In that case,
                       'trimOverlap=True' can be used to ensure that a given
                       interval is never seen more than once.

        numericGroups: Whether to return group labels or simply the numeric
                       index. The latter is more useful when these are passed to
                       a function whose output will be sorted according to group.

        includeStrand: Whether to include the strand in the output. The default is False

        >>> from deeptoolsintervals import parse
        >>> from os.path import dirname, basename
        >>> gtf = parse.GTF(["{0}/test/GRCh38.84.bed6".format(dirname(parse.__file__)), "{0}/test/GRCh38.84.bed2".format(dirname(parse.__file__))], keepExons=True)
        >>> overlaps = gtf.findOverlaps("1", 0, 3000000)
        >>> labels = dict()
        >>> for o in overlaps:
        ...     if basename(o[3]) not in labels:
        ...         labels[basename(o[3])] = 0
        ...     labels[basename(o[3])] += 1
        >>> assert(labels['GRCh38.84.bed2'] == 1)
        >>> assert(labels['GRCh38.84.bed6'] == 15)
        >>> assert(labels['group2'] == 9)
        >>> assert(labels['group 3'] == 7)
        >>> assert(labels['group 1'] == 6)
        >>> assert(labels['group 1_r1'] == 4)
        """
        chrom = self.mungeChromosome(chrom, append=False)
        if not chrom:
            return None

        # Ensure that this is a tree and has entries, otherwise
        if self.tree.countEntries() == 0:
            return None

        if not self.tree.isTree():
            raise RuntimeError('The GTFtree is actually a vine! There must have been an error during creation (this shouldn\'t happen)...')

        # Convert the strand to a number
        if strand == '+':
            strand = 1
        elif strand == '-':
            strand = 2
        else:
            strand = 0

        overlaps = self.tree.findOverlaps(chrom, start, end, strand, matchType, strandType, "transcript_id", includeStrand)
        if overlaps is None:
            return None

        for i, o in enumerate(overlaps):
            if o[2] not in self.exons or len(self.exons[o[2]]) == 0:
                exons = [(o[0], o[1])]
            else:
                exons = sorted(self.exons[o[2]])

            if numericGroups:
                overlaps[i] = (o[0], o[1], o[2], o[3], exons)
            else:
                overlaps[i] = (o[0], o[1], o[2], self.labels[o[3]], exons)

            if includeStrand:
                overlaps[i] = overlaps[i] + (str(o[-2].decode("ascii")),)

            # Add the score
            overlaps[i] = overlaps[i] + (o[-1],)

        # Ensure that the intervals are sorted by their 5'-most bound. This enables trimming
        overlaps = sorted(overlaps)

        if trimOverlap:
            while True:
                if len(overlaps) > 0:
                    if overlaps[0][0] < start:
                        del overlaps[0]
                    else:
                        break
                else:
                    overlaps = []
                    break

        return overlaps
