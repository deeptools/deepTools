#!/usr/bin/env python
import xmlrpclib
import time
import numpy as np

class deepBlue(object):
    def __init__(self, sample, url="http://deepblue.mpi-inf.mpg.de/xmlrpc", userKey="anonymous_key"):
        """
        Connect to the requested deepblue server with the given user key and request the specifed sample from it.
        """
        self.sample = sample
        self.url = url
        self.userKey = userKey
        self.server = xmlrpclib.Server(url, allow_none=True)
        self.info = None
        self.sampleID = None
        self.genome = None
        self.chroms = None

        # Set self.sampleID
        id = self.sample2id()
        if not id:
            raise RuntimeError("The requested sample ({}) is not available".format(sample))

        # Set self.info
        (status, resp) = server.info(id, userKey)
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

    def sample2id(self):
        """
        Given a sample name, return its ID (or None on error). On error, this raises a runtime exception.
        
        self.sampleID is then the internal sample ID (e.g., s9971)
        """
        (status, resp) = self.server.search(self.sample, "samples", self.userKey)
        if status != "okay":
            raise RuntimeError("Received an error ({}) while searching for '{}'".format(status, sample))
        if len(resp) == 0:
            raise RuntimeError("'{}' does not exist".format(self.sample))
        if len(resp) > 1:
            raise RuntimeError("'{}' is an ambiguous name".format(self.sample))
        self.sampleID = resp[0][0]
        return self.sampleID
    
    def getGenome(self):
        """
        Determines and sets the genome assigned to a given sample. On error, this raises a runtime exception.
        
        self.genome is then the internal genome ID.
        """
        (status, resp) = self.server.search(self.sample, "genomes", self.userKey)
        if status != "okay":
            raise RuntimeError("Searching for the genome for '{}' caused an error: {}".format(sample, status))
        if len(resp) == 0:
            raise RuntimeError("No annotated genome for '{}'".format(sample))
        if len(resp) > 1:
            raise RuntimeError("Multiple possible genomes for '{}'!".format(sample))

        self.genome = resp[0][0]
        return self.genome
    
    def getChroms(self):
        """
        Determines and sets the chromosome names/sizes for a given sample. On error, this raises a runtime exception.
        
        self.chroms is then a dictionary of chromosome:length pairs
        """
        (status, resp) = self.server.chromosomes(self.genome, self.userKey)
        if status != "okay":
            raise RuntimeError("Received an error while fetching chromosome information for '{}': {}".format(sample, resp))
        self.chroms = {k: v for k, v in resp}
        return resp

    def getValuesInRegion(self, chrom, start=0, end=0):
        """
        Return a numpy array with values for each base in the requested region.
        If end==0 then the region from start to the end of the chromosome is returned.
        """
        if chrom not in self.chroms.keys():
            raise RuntimeError("'{}' is not a valid chromosome.".format(chrom))
    
        if end == 0:
            end = self.chroms[chrom]

        if start >= end:
            raise RuntimeError("The start position MUST be less then the end position ({} and {})".format(start, end))

        # Get the experiment information
        (status, queryID) = self.server.select_experiments(self.sample, chrom, start, end, self.userKey)
        if status != "okay":
            raise RuntimeError("Received the following error while fetching values in the range {}:{}-{} in file '{}': {}".format(chrom, start, end, self.sample, queryID))
        if not queryID:
            raise RuntimeError("Somehow, we received None as a query ID (range {}:{}-{} in file '{}')".format(chrom, start, end, self.sample))

        # Query the regions
        (status, reqID) = self.server.get_regions(queryID, "START,END,VALUE", self.userKey)
        if status != "okay":
            raise RuntimeError("Received the following error while fetching regions in the range {}:{}-{} in file '{}': {}".format(chrom, start, end, self.sample, reqID))

        # Wait for the server to process the data
        (status, info) = self.server.info(reqID, self.userKey)
        request_status = info[0]["state"]
        while request_status != "done" and request_status != "failed":
            time.sleep(1)
            (status, info) = self.server.info(reqID, self.userKey)
            request_status = info[0]["state"]

        # Get the actual data
        (status, resp) = self.server.get_request_data(reqID, self.userKey)
        if status != "okay":
            raise RuntimeError("Received the following error while fetching data in the range {}:{}-{} in file '{}': {}".format(chrom, start, end, self.sample, resp))

        # Generate the output
        o = np.empty(end-start)
        o[:] = np.nan
        for intervals in resp.split("\n"):
            interval = intervals.split("\t")
            s = int(interval[0]) - start
            e = s + int(interval[1]) - int(interval[0])
            o[s:e] = float(interval[2])

        return o
