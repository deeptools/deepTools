#!/usr/bin/env python

from deeptoolsintervals import tree
from deeptoolsintervals.parse import GTF, openPossiblyCompressed
import sys
from os.path import basename
import csv


class Enrichment(GTF):
    """
    This is like the GTF object, but has no groups or exons (but a "features" list). BED files are given a 'peaks' feature and GTF files use column 3.
    """

    def parseBEDcore(self, line, ncols, feature):
        strand = 3
        cols = line.split("\t")

        if int(cols[1]) < 0:
            cols[1] = 0

        if int(cols[1]) >= int(cols[2]):
            sys.stderr.write("Warning: {0}:{1}-{2} is an invalid BED interval! Ignoring it.\n".format(cols[0], cols[1], cols[2]))
            return

        # BED6/BED12: set name and strand
        score = '.'
        if ncols > 3:
            if cols[5] == '+':
                strand = 0
            elif cols[5] == '-':
                strand = 1
            score = cols[4]

        if ncols != 12 or self.keepExons is False:
            self.tree.addEnrichmentEntry(self.mungeChromosome(cols[0]), int(cols[1]), int(cols[2]), strand, score, feature)
        else:
            starts = cols[11].strip(",").split(",")
            widths = cols[10].strip(",").split(",")
            starts = [int(x) + int(cols[1]) for x in starts]
            ends = [x + int(y) for x, y in zip(starts, widths)]
            for x, y in zip(starts, ends):
                self.tree.addEnrichmentEntry(self.mungeChromosome(cols[0]), x, y, strand, score, feature)

    def parseBED(self, fp, line, ncols=3, feature='peaks', labelColumn=None):
        """
        parse a BED file. The default feature label is 'peaks'

        fp:    A python file pointer
        line:  The first line
        ncols: The number of columns to care about
        feature: The feature label
        labelColumn: If this isn't None, it overrides the 'feature' option

        >>> from deeptoolsintervals import enrichment
        >>> from os.path import dirname
        >>> gtf = enrichment.Enrichment("{0}/test/GRCh38.84.bed".format(dirname(enrichment.__file__)), keepExons=True)
        >>> o = gtf.findOverlaps("1", [(1, 3000000)])
        >>> assert(o == frozenset(['GRCh38.84.bed']))
        >>> o = gtf.findOverlaps("chr1", [(1, 3000000)])
        >>> assert(o == frozenset(['GRCh38.84.bed']))
        """

        # Handle the first line
        if labelColumn is not None:
            cols = line.split("\t")
            feature = cols.pop(labelColumn)
            line = "\t".join(cols)
        self.parseBEDcore(line, ncols, feature)
        if feature not in self.features:
            self.features.append(feature)

        # iterate over the remaining lines
        for line in fp:
            if not isinstance(line, str):
                line = line.decode('ascii')
            line = line.strip()
            if len(line) == 0:
                # Apparently this happens, some people seem to like trying to break things
                continue

            if line.startswith("#"):
                continue
            else:
                if labelColumn is not None:
                    cols = line.split("\t")
                    feature = cols.pop(labelColumn)
                    line = "\t".join(cols)
                self.parseBEDcore(line, ncols, feature)

            if feature not in self.features:
                self.features.append(feature)

    def parseGTF(self, fp, line):
        """
        >>> from deeptoolsintervals import enrichment
        >>> from os.path import dirname
        >>> gtf = enrichment.Enrichment("{0}/test/GRCh38.84.gtf.gz".format(dirname(enrichment.__file__)), keepExons=True)
        >>> o = gtf.findOverlaps("1", [(0, 2000000)])
        >>> assert(o == frozenset(['start_codon', 'exon', 'stop_codon', 'CDS', 'gene', 'transcript', 'group 1', 'group 2']))
        """

        # Handle the first line
        cols = line.split("\t")

        strand = 3
        if cols[6] == '+':
            strand = 0
        elif cols[6] == '-':
            strand = 1

        feature = cols[2]
        if "deepTools_group" in cols[8]:
            s = next(csv.reader([cols[8]], delimiter=' '))
            if s[-1] != "deepTools_group":
                feature = s[s.index("deepTools_group") + 1].rstrip(";")

        self.tree.addEnrichmentEntry(self.mungeChromosome(cols[0]), int(cols[3]) - 1, int(cols[4]), strand, cols[5], feature)
        if feature not in self.features:
            self.features.append(feature)

        # Handle the remaining lines
        for line in fp:
            if not isinstance(line, str):
                line = line.decode('ascii')
            if not line.startswith('#'):
                cols = line.split("\t")
                if len(cols) == 0:
                    continue

                strand = 3
                if cols[6] == '+':
                    strand = 0
                elif cols[6] == '-':
                    strand = 1

                feature = cols[2]
                if "deepTools_group" in cols[8]:
                    s = next(csv.reader([cols[8]], delimiter=" "))
                    if s[-1] != "deepTools_group":
                        feature = s[s.index("deepTools_group") + 1].rstrip(";")

                self.tree.addEnrichmentEntry(self.mungeChromosome(cols[0]), int(cols[3]) - 1, int(cols[4]), strand, cols[5], feature)
                if feature not in self.features:
                    self.features.append(feature)

    def __init__(self, fnames, keepExons=False, labels=None, verbose=False):
        """
        Driver function to actually parse files. The steps are as follows:

        1) skip to the first non-comment line
        2) Infer the type from that
        3) Call a type-specific processing function accordingly
          * These call the underlying C code for storage
          * These handle chromsome name conversions (python-level)

        Required inputs are as follows:

        fnames:	A list of (possibly compressed with gzip or bzip2) GTF or BED files.

        Optional input is:

        keepExons:    For BED12 files, exons are ignored by default.
        labels:       Override the feature labels supplied in the file(s).
                      Note that this might instead be replaced later in the .features attribute.
        verbose:      Whether to print warnings (default: False)
        """
        self.fname = []
        self.filename = ""
        self.chroms = []
        self.features = []
        self.tree = tree.initTree()
        self.keepExons = keepExons
        self.verbose = verbose

        if not isinstance(fnames, list):
            fnames = [fnames]

        # Load the files
        for labelIdx, fname in enumerate(fnames):
            self.filename = fname
            fp = openPossiblyCompressed(fname)
            line, labelColumn = self.firstNonComment(fp)
            if line is None:
                # This will only ever happen if a file is empty or just has a header/comment
                continue
            line = line.strip()

            ftype = self.inferType(fp, line, labelColumn)

            if ftype != 'GTF' and labels is not None:
                assert(len(labels) > labelIdx)
                bname = labels[labelIdx]
            else:
                bname = basename(fname)
            if ftype == 'GTF':
                self.parseGTF(fp, line)
            elif ftype == 'BED3':
                self.parseBED(fp, line, 3, feature=bname, labelColumn=labelColumn)
            elif ftype == 'BED6':
                self.parseBED(fp, line, 6, feature=bname, labelColumn=labelColumn)
            else:
                self.parseBED(fp, line, 12, feature=bname, labelColumn=labelColumn)
            fp.close()

        # Sanity check
        if self.tree.countEntries() == 0:
            raise RuntimeError("None of the input BED/GTF files had valid regions")

        if len(self.features) == 0:
            raise RuntimeError("There were no valid feature labels!")

        # vine -> tree
        self.tree.finish()

    # findOverlaps()
    def findOverlaps(self, chrom, blocks, strand=".", matchType=0, strandType=0):
        """
        Given a chromosome and start/end coordinates with an optional strand,
        return a frozenset of the overlap features.

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
        """
        chrom = self.mungeChromosome(chrom, append=False)
        if not chrom:
            return None

        # Ensure that this is a tree and has entries
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

        oset = frozenset()
        for block in blocks:
            overlaps = self.tree.findOverlappingFeatures(chrom, int(block[0]), int(block[1]), strand, matchType, strandType)
            if overlaps is not None:
                oset = oset.union(frozenset(overlaps))

        return oset
