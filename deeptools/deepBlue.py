#!/usr/bin/env python
try:
    # python 2
    import xmlrpclib
except:
    # python 3
    import xmlrpc.client as xmlrpclib
import time
import tempfile
import os.path
import sys
import pyBigWig
from deeptools.utilities import mungeChromosome
from deeptoolsintervals import GTF
import datetime


def isDeepBlue(fname):
    """
    Returns true if the file ends in .wig, .wiggle, or .bedgraph, since these indicate a file on the deepBlue server
    """
    if fname.endswith(".wig"):
        return True
    if fname.endswith(".wiggle"):
        return True
    if fname.endswith(".bedgraph"):
        return True
    if fname.startswith("http") or fname.startswith("ftp"):
        return False
    # For ENCODE samples, the "Name" is just the ENCODE sample ID, so as a fallback check for files that aren't there.
    if not os.path.exists(fname):
        return True
    return False


def mergeRegions(regions):
    """
    Given a list of [(chrom, start, end), ...], merge all overlapping regions

    This returns a dict, where values are sorted lists of [start, end].
    """
    bar = sorted(regions)
    out = dict()
    last = [None, None, None]
    for reg in bar:
        if reg[0] == last[0] and reg[1] <= last[2]:
            if reg[2] > last[2]:
                last[2] = reg[2]
            continue
        else:
            if last[0]:
                if last[0] not in out:
                    out[last[0]] = list()
                out[last[0]].append([last[1], last[2]])
            last = [reg[0], reg[1], reg[2]]
    if last[0] not in out:
        out[last[0]] = list()
    out[last[0]].append([last[1], last[2]])
    return out


def makeTiles(db, args):
    """
    Given a deepBlue object, return a list of regions that will be queried
    """
    out = []
    for (k, v) in db.chromsTuple:
        start = 0
        while start <= v:
            end = start + args.binSize
            if end > v:
                end = v
            out.append([k, start, end])
            start += end + args.distanceBetweenBins
    return out


def makeChromTiles(db):
    """
    Make a region for each chromosome
    """
    out = []
    for (k, v) in db.chromsTuple:
        out.append([k, 0, v])
    return out


def makeRegions(BED, args):
    """
    Given a list of BED/GTF files, make a list of regions.
    These are vaguely extended as appropriate. For simplicity, the maximum of --beforeRegionStartLength
    and --afterRegionStartLength are tacked on to each end and transcripts are used for GTF files.
    """
    itree = GTF(BED, transcriptID=args.transcriptID, transcript_id_designator=args.transcript_id_designator)
    o = []
    extend = 0
    # The before/after stuff is specific to computeMatrix
    if "beforeRegionStartLength" in args:
        extend = max(args.beforeRegionStartLength, args.afterRegionStartLength)
    for chrom in itree.chroms:
        regs = itree.findOverlaps(chrom, 0, 4294967295)  # bigWig files use 32 bit coordinates
        for reg in regs:
            o.append([chrom, max(0, reg[0] - extend), reg[1] + extend])
    del itree
    return o


def preloadWrapper(foo):
    """
    This is a wrapper around the preload function for multiprocessing
    """
    args = foo[2]
    regs = foo[3]
    res = deepBlue(foo[0], url=args.deepBlueURL, userKey=args.userKey)
    return res.preload(regs, tmpDir=args.deepBlueTempDir)


class deepBlue(object):
    def __init__(self, sample, url="http://deepblue.mpi-inf.mpg.de/xmlrpc", userKey="anonymous_key"):
        """
        Connect to the requested deepblue server with the given user key and request the specifed sample from it.

        >>> sample = "S002R5H1.ERX300721.H3K4me3.bwa.GRCh38.20150528.bedgraph"
        >>> db = deepBlue(sample) # doctest: +SKIP
        >>> assert(db.chroms("chr1") == 248956422) # doctest: +SKIP
        """
        self.sample = sample
        self.url = url
        self.userKey = userKey
        self.server = xmlrpclib.Server(url, allow_none=True)
        self.info = None
        self.experimentID = None
        self.genome = None
        self.chromsDict = None
        self.chromsTuple = None

        # Set self.experimentID
        experimentID = self.getEID()
        if not experimentID:
            raise RuntimeError("The requested sample({}) has no associated experiment! If you did not intend to use samples on deepBlue, then it appears either you misspelled a file name or (if you're using BAM files for input) one of your BAM files is lacking a valid index.".format(sample))

        # Set self.info
        (status, resp) = self.server.info(self.experimentID, userKey)
        if status != "okay":
            raise RuntimeError("Received the following error while fetching information about '{}': {}".format(resp, sample))
        self.info = resp[0]

        # Set self.genome
        genome = self.getGenome()
        if not genome:
            raise RuntimeError("Unable to determine an appropriate genome for '{}'".format(sample))

        # Set self.chroms
        chroms = self.getChroms()
        if not chroms:
            raise RuntimeError("Unable to determine chromosome names/sizes for '{}'".format(sample))

    def getEID(self):
        """
        Given a sample name, return its associated experiment ID (or None on error).

        self.experimentID is then the internal ID (e.g., e52525)
        """
        (status, resps) = self.server.search(self.sample, "experiments", self.userKey)
        if status != "okay":
            raise RuntimeError("Received an error ({}) while searching for the experiment associated with '{}'".format(resps, self.sample))
        for resp in resps:
            if resp[1] == self.sample:
                self.experimentID = resp[0]
                return resp[0]
        return None

    def getGenome(self):
        """
        Determines and sets the genome assigned to a given sample. On error, this raises a runtime exception.

        self.genome is then the internal genome ID.
        """
        if "genome" in self.info.keys():
            self.genome = self.info["genome"]
        return self.genome

    def getChroms(self):
        """
        Determines and sets the chromosome names/sizes for a given sample. On error, this raises a runtime exception.

        self.chroms is then a dictionary of chromosome:length pairs
        """
        (status, resp) = self.server.chromosomes(self.genome, self.userKey)
        if status != "okay":
            raise RuntimeError("Received an error while fetching chromosome information for '{}': {}".format(self.sample, resp))
        self.chromsDict = {k: v for k, v in resp}
        self.chromsTuple = [(k, v) for k, v in resp]
        return resp

    def chroms(self, chrom=None):
        """
        Like the chroms() function in pyBigWig, returns either chromsDict (chrom is None) or the length of a given chromosome
        """
        if chrom is None:
            return self.chromsDict
        elif chrom in self.chromsDict:
            return self.chromsDict[chrom]
        return None

    def close(self):
        pass

    def preload(self, regions, tmpDir=None):
        """
        Given a sample and a set of regions, write a bigWig file containing the underlying signal.

        This function returns the file name, which needs to be deleted by the calling function at some point.

        This sends queries one chromosome at a time, due to memory limits on deepBlue
        """
        startTime = datetime.datetime.now()
        regions2 = mergeRegions(regions)

        # Make a temporary file
        f = tempfile.NamedTemporaryFile(delete=False, dir=tmpDir)
        fname = f.name
        f.close()

        # Start with the bigWig file
        bw = pyBigWig.open(fname, "w")
        bw.addHeader(self.chromsTuple, maxZooms=0)  # This won't work in IGV!

        # Make a string out of everything in a resonable order
        for k, v in self.chromsTuple:
            # Munge chromosome names as appropriate
            chrom = mungeChromosome(k, regions2.keys())
            if not chrom:
                continue
            if chrom not in regions2 or len(regions2) == 0:
                continue
            regionsStr = "\n".join(["{}\t{}\t{}".format(k, reg[0], reg[1]) for reg in regions2[chrom]])
            regionsStr += "\n"

            # Send the regions
            (status, regionsID) = self.server.input_regions(self.genome, regionsStr, self.userKey)
            if status != "okay":
                raise RuntimeError("Received the following error while sending regions for '{}': {}".format(regionsID, self.sample))

            # Get the experiment information
            (status, queryID) = self.server.select_experiments(self.sample, k, None, None, self.userKey)
            if status != "okay":
                raise RuntimeError("Received the following error while running select_experiments on file '{}': {}".format(self.sample, queryID))
            if not queryID:
                raise RuntimeError("Somehow, we received None as a query ID (file '{}')".format(self.sample))

            # Intersect
            (status, intersectID) = self.server.intersection(queryID, regionsID, self.userKey)
            if status != "okay":
                raise RuntimeError("Received the following error while running intersection on file '{}': {}".format(self.sample, intersectID))
            if not intersectID:
                raise RuntimeError("Somehow, we received None as an intersect ID (file '{}')".format(self.sample))

            # Query the regions
            (status, reqID) = self.server.get_regions(intersectID, "START,END,VALUE", self.userKey)
            if status != "okay":
                raise RuntimeError("Received the following error while fetching regions in file '{}': {}".format(self.sample, reqID))

            # Wait for the server to process the data
            (status, info) = self.server.info(reqID, self.userKey)
            request_status = info[0]["state"]
            while request_status != "done" and request_status != "failed":
                time.sleep(0.1)
                (status, info) = self.server.info(reqID, self.userKey)
                request_status = info[0]["state"]

            # Get the actual data
            (status, resp) = self.server.get_request_data(reqID, self.userKey)
            if status != "okay":
                raise RuntimeError("Received the following error while fetching data in file '{}': {}".format(self.sample, resp))

            for intervals in resp.split("\n"):
                interval = intervals.split("\t")
                if interval[0] == '':
                    continue
                bw.addEntries([k], [int(interval[0]) - 1], ends=[int(interval[1]) - 1], values=[float(interval[2])])
        bw.close()
        sys.stderr.write("{} done (took {})\n".format(self.sample, datetime.datetime.now() - startTime))
        sys.stderr.flush()

        return fname
