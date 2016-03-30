#!/usr/bin/env python

from deeptoolsintervals import tree
import re
import sys
import gzip
import bz2
import os.path


def getNext(fp):
    """
    Sometimes we need to decode, sometimes not
    """
    line = fp.readline()
    if isinstance(line, str):
        return line
    return line.decode('ascii')


def seemsLikeGTF(cols, regex):
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
        assert(regex.match(cols[8]) is not None)
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
        name = name + "_r" + str(i)
        if name not in labels:
            return name


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
    elif first3 == b"\x42\x5a\x68":
        return bz2.BZ2File(fname, "rb")
    else:
        return open(fname)


class GTF(object):
    """
    A class to hold an interval tree and its associated functions

    >>> from deeptoolsintervals import parse
    >>> from os.path import dirname
    >>> gtf = parse.GTF("{0}/test/GRCh38.84.gtf.gz".format(dirname(parse.__file__)), keepExons=True)
    >>> gtf.findOverlaps("1", 1, 20000)
    [(11868, 14409, 'ENST00000456328', 'group 1', [(11868, 12227), (12612, 12721), (13220, 14409)]), (12009, 13670, 'ENST00000450305', 'group 1', [(12009, 12057), (12178, 12227), (12612, 12697), (12974, 13052), (13220, 13374), (13452, 13670)]), (14403, 29570, 'ENST00000488147', 'group 1', [(14403, 14501), (15004, 15038), (15795, 15947), (16606, 16765), (16857, 17055), (17232, 17368), (17605, 17742), (17914, 18061), (18267, 18366), (24737, 24891), (29533, 29570)]), (17368, 17436, 'ENST00000619216', 'group 2', [(17368, 17436)])]
    >>> gtf = parse.GTF("{0}/test/GRCh38.84.gtf.gz".format(dirname(parse.__file__)))
    >>> gtf.findOverlaps("1", 1, 20000)
    [(11868, 14409, 'ENST00000456328', 'group 1', [(11868, 14409)]), (12009, 13670, 'ENST00000450305', 'group 1', [(12009, 13670)]), (14403, 29570, 'ENST00000488147', 'group 1', [(14403, 29570)]), (17368, 17436, 'ENST00000619216', 'group 2', [(17368, 17436)])]
    >>> gtf.findOverlaps("1", 12000, 20000, trimOverlap=True)
    [(12009, 13670, 'ENST00000450305', 'group 1', [(12009, 13670)]), (14403, 29570, 'ENST00000488147', 'group 1', [(14403, 29570)]), (17368, 17436, 'ENST00000619216', 'group 2', [(17368, 17436)])]
    >>> gtf.findOverlaps("1", 1, 20000, numericGroups=True, includeStrand=True)
    [(11868, 14409, 'ENST00000456328', 0, [(11868, 14409)], '+'), (12009, 13670, 'ENST00000450305', 0, [(12009, 13670)], '+'), (14403, 29570, 'ENST00000488147', 0, [(14403, 29570)], '-'), (17368, 17436, 'ENST00000619216', 1, [(17368, 17436)], '-')]
    """

    def firstNonComment(self, fp):
        """
        Skip lines at the beginning of a file starting with #, browser, or track.
        """
        line = getNext(fp)
        try:
            while line.startswith("#") or line.startswith('track') or line.startswith('browser'):
                line = getNext(fp)
        except:
            sys.stderr.write("Warning, {0} was empty\n".format(self.filename))
            return None
        return line

    def inferType(self, fp, line):
        """
        Attempt to infer a file type from a single line. This is largely based on the number of columns plus a regex looking for "gene_id".
        """
        cols = line.split("\t")
        if len(cols) < 3:
            raise RuntimeError('{0} does not seem to be a recognized file type!'.format(self.filename))
        elif len(cols) == 3:
            return 'BED3'
        elif len(cols) < 6:
            sys.stderr.write("Warning, {0} has an abnormal number of fields. Assuming BED3 format.\n".format(self.filename))
            return 'BED3'
        elif len(cols) == 6:
            return 'BED6'
        elif len(cols) == 9 and seemsLikeGTF(cols, self.gene_id_regex):
            return 'GTF'
        elif len(cols) == 12:
            return 'BED12'
        elif len(cols) < 12:
            sys.stderr.write("Warning, {0} has an abnormal format. Assuming BED6 format.\n".format(self.filename))
            return 'BED6'
        else:
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
        else:
            return None

        return chrom

    def parseBEDcore(self, line, ncols):
        """
        Returns True if the entry was added, otherwise False

        >>> from deeptoolsintervals import parse
        >>> from os.path import dirname
        >>> gtf = parse.GTF("{0}/test/GRCh38.84.bed12.bz2".format(dirname(parse.__file__)), keepExons=True, labels=["foo"])
        >>> gtf.findOverlaps("1", 1, 20000)
        [(11868, 14409, 'ENST00000456328.2', 'foo', [(11868, 12227), (12612, 12721), (13220, 14409)]), (12009, 13670, 'ENST00000450305.2', 'foo', [(12009, 12057), (12178, 12227), (12612, 12697), (12974, 13052), (13220, 13374), (13452, 13670)]), (14403, 29570, 'ENST00000488147.1', 'foo', [(14403, 14501), (15004, 15038), (15795, 15947), (16606, 16765), (16857, 17055), (17232, 17368), (17605, 17742), (17914, 18061), (18267, 18366), (24737, 24891), (29533, 29570)]), (17368, 17436, 'ENST00000619216.1', 'foo', [(17368, 17436)])]
        """

        strand = 3
        cols = line.split("\t")
        name = "{0}:{1}-{2}".format(cols[0], cols[1], cols[2])

        if int(cols[1]) >= int(cols[2]) or int(cols[1]) < 0:
            sys.stderr.write("Warning: {0}:{1}-{2} is an invalid BED interval! Ignoring it.\n".format(cols[0], cols[1], cols[2]))
            return

        # BED6/BED12: set name and strand
        if ncols > 3:
            name = cols[3]
            if cols[5] == '+':
                strand = 0
            elif cols[5] == '-':
                strand = 1

        # Ensure that the name is unique
        if name in self.exons.keys():
            sys.stderr.write("Skipping {0}, an entry by this name already exists!\n".format(name))
            return False
        else:
            self.tree.addEntry(self.mungeChromosome(cols[0]), int(cols[1]), int(cols[2]), name, strand, self.labelIdx)
            if ncols != 12 or self.keepExons is False:
                self.exons[name] = [(int(cols[1]), int(cols[2]))]
            else:
                assert(len(cols) == 12)
                self.exons[name] = parseExonBounds(int(cols[1]), int(cols[2]), int(cols[9]), cols[10], cols[11])
        return True

    def parseBED(self, fp, line, ncols=3):
        """
        parse a BED file. The default group label is the file name.

        fp:    A python file pointer
        line:  The first line
        ncols: The number of columns to care about


        >>> from deeptoolsintervals import parse
        >>> from os.path import dirname
        >>> gtf = parse.GTF("{0}/test/GRCh38.84.bed6".format(dirname(parse.__file__)), keepExons=True)
        >>> gtf.findOverlaps("1", 1, 20000)
        [(11868, 14409, 'ENST00000456328.2', 'group 1', [(11868, 14409)]), (12009, 13670, 'ENST00000450305.2', 'group 1', [(12009, 13670)]), (14403, 29570, 'ENST00000488147.1', 'group 1', [(14403, 29570)]), (17368, 17436, 'ENST00000619216.1', 'group 1', [(17368, 17436)])]
        >>> gtf = parse.GTF("{0}/test/GRCh38.84.bed".format(dirname(parse.__file__)), keepExons=True, labels=["foo", "bar", "quux", "sniggly"])
        >>> gtf.findOverlaps("1", 1, 20000)
        [(11868, 14409, '1:11868-14409', 'foo', [(11868, 14409)]), (12009, 13670, '1:12009-13670', 'foo', [(12009, 13670)]), (14403, 29570, '1:14403-29570', 'foo', [(14403, 29570)]), (17368, 17436, '1:17368-17436', 'foo', [(17368, 17436)])]
        """
        groupLabelsFound = 0
        groupEntries = 0

        # Handle the first line
        if self.parseBEDcore(line, ncols):
            groupEntries = 1

        # iterate over the remaining lines
        for line in fp:
            if not isinstance(line, str):
                line = line.decode('ascii')
            line = line.strip()
            if len(line) == 0:
                # Apparently this happens, some people seem to like trying to break things
                continue

            if line.startswith("#"):
                # If there was a previous group AND it had no entries then remove it
                if groupLabelsFound > 0:
                    if groupEntries == 0:
                        sys.stderr.write("Warning, the '{0}' group had no valid entries! Removing it.\n".format(self.labels[self.labelIdx]))
                        del self.labels[-1]
                        groupLabelsFound -= 1
                        self.labelIdx -= 1

                label = line[1:].strip()
                if len(label):
                    # Gaard against duplicate group labels
                    self.labels.append(findRandomLabel(self.labels, label))
                else:
                    # I'm sure someone will try an empty label...
                    self.labels.append(findRandomLabel(self.labels, os.path.basename(self.filename)))
                self.labelIdx += 1
                groupLabelsFound += 1
                groupEntries = 0
            else:
                if self.parseBEDcore(line, ncols):
                    groupEntries += 1

        if groupEntries > 0:
            if self.defaultGroup is not None:
                self.labels.append(findRandomLabel(self.labels, self.defaultGroup))
            else:
                self.labels.append(findRandomLabel(self.labels, os.path.basename(self.filename)))
            self.labelIdx += 1

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

        m = self.deepTools_group_regex.search(cols[8])
        if m:
            label = m.groups()[0]
        elif self.defaultGroup is not None:
            label = self.defaultGroup

        m = self.transcript_id_regex.search(cols[8])
        if not m:
            sys.stderr.write("Warning: {0} is malformed!\n".format("\t".join(cols)))
            return

        name = m.groups()[0]
        if name in self.exons.keys():
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

        # Get the label index
        if label not in self.labels:
            self.labels.append(label)
        self.labelIdx = self.labels.index(label)

        chrom = self.mungeChromosome(cols[0])
        self.tree.addEntry(chrom, int(cols[3]) - 1, int(cols[4]), name, strand, self.labelIdx)

        # Exon bounds placeholder
        self.exons[name] = []

    def parseGTFexon(self, cols):
        """
        Parse an exon entry and add it to the transcript hash
        """
        if int(cols[3]) - 1 < 0:
            sys.stderr.write("Warning: Invalid start in '{0}', skipping\n".format("\t".join(cols)))
            return

        m = self.transcript_id_regex.search(cols[8])
        if not m:
            sys.stderr.write("Warning: {0} is malformed!\n".format("\t".join(cols)))
            return

        name = m.groups()[0]
        if name in self.transcriptIDduplicated:
            return
        if name not in self.exons.keys():
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
        ...     if basename(o[3]) not in labels.keys():
        ...         labels[basename(o[3])] = 0
        ...     labels[basename(o[3])] += 1
        >>> assert(labels['GRCh38.84.gtf.gz'] == 17)
        >>> assert(labels['GRCh38.84.2.gtf.gz'] == 6)
        >>> assert(labels['group 1'] == 5)
        >>> assert(labels['group 2'] == 3)
        """
        file_label = findRandomLabel(self.labels, os.path.basename(self.filename))

        # Handle the first line
        cols = line.split("\t")
        if cols[2].lower() == self.transcriptID:
            self.parseGTFtranscript(cols, file_label)
        elif cols[2].lower() == self.exonID:
            self.parseGTFexon(cols)

        # Handle the remaining lines
        for line in fp:
            if not isinstance(line, str):
                line = line.decode('ascii')
            if not line.startswith('#'):
                cols = line.split("\t")
                if len(cols) == 0:
                    continue

                if cols[2].lower() == self.transcriptID:
                    self.parseGTFtranscript(cols, file_label)
                elif cols[2].lower() == self.exonID and self.keepExons is True:
                    self.parseGTFexon(cols)

        # Reset self.labelIdx
        self.labelIdx = len(self.labels) - 1

    def __init__(self, fnames, exonID="exon", transcriptID="transcript", keepExons=False, labels=[], transcript_id_designator="transcript_id", defaultGroup=None):
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
        transcript_id_designator: For gtf files, this the key used in a regex.
                      If one set transcriptID to 'gene', then
                      transcript_id_designator would need to be changed to
                      'gene_id' or 'gene_name' to extract the gene ID/name from
                      the attributes.
        defaultGroup: The default group name. If None, the file name is used.
        """
        self.fname = []
        self.filename = ""
        self.chroms = []
        self.exons = {}
        self.labels = []
        self.transcriptIDduplicated = []
        self.tree = tree.initTree()
        self.labelIdx = 0
        self.gene_id_regex = re.compile('(?:gene_id (?:\"([ \w\d"\-]+)\"|([ \w\d"\-]+))[;|\r|\n])')
        self.transcript_id_regex = re.compile('(?:{0} (?:\"([ \w\d"\-]+)\"|([ \w\d"\-]+))[;|\r|\n])'.format(transcript_id_designator))
        self.deepTools_group_regex = re.compile('(?:deepTools_group (?:\"([ \w\d"\-]+)\"|([ \w\d"\-]+))[;|\r|\n])')
        self.exonID = exonID
        self.transcriptID = transcriptID
        self.keepExons = keepExons
        self.defaultGroup = defaultGroup

        if labels != []:
            self.already_input_labels = True

        if not isinstance(fnames, list):
            fnames = [fnames]

        # Load the files
        for fname in fnames:
            self.filename = fname
            fp = openPossiblyCompressed(fname)
            line = self.firstNonComment(fp).strip()
            if line is None:
                # This will only ever happen if a file is empty or just has a header/comment
                continue

            ftype = self.inferType(fp, line)
            if ftype == 'GTF':
                self.parseGTF(fp, line)
            elif ftype == 'BED3':
                self.parseBED(fp, line, 3)
            elif ftype == 'BED6':
                self.parseBED(fp, line, 6)
            else:
                self.parseBED(fp, line, 12)
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
        ...     if basename(o[3]) not in labels.keys():
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
        if not overlaps:
            return None

        for i, o in enumerate(overlaps):
            if o[2] not in self.exons.keys() or len(self.exons[o[2]]) == 0:
                exons = [(o[0], o[1])]
            else:
                exons = sorted(self.exons[o[2]])

            if numericGroups:
                overlaps[i] = (o[0], o[1], o[2], o[3], exons)
            else:
                overlaps[i] = (o[0], o[1], o[2], self.labels[o[3]], exons)

            if includeStrand:
                overlaps[i] = overlaps[i] + (str(o[-1].decode("ascii")),)

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
                    overlaps = None
                    break

        return overlaps
