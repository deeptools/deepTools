#from unittest import TestCase

import deeptools.countReadsPerBin as cr
import numpy as np
import numpy.testing as nt

__author__ = 'Fidel'


class TestCountReadsPerBin(object):
    def setUp(self):
        """
        The distribution of reads between the two bam files is as follows.

        They cover 200 bp::

              0                              100                           200
              |------------------------------------------------------------|
            A                                ===============
                                                            ===============


            B                 ===============               ===============
                                             ===============
                                                            ===============
        """
        self.root = "./test/test_data/"
        self.bamFile1  = self.root + "testA.bam"
        self.bamFile2  = self.root + "testB.bam"
        self.bamFile_PE  = self.root + "test_paired2.bam"
        self.chrom = '3R'
        step_size = 50
        bin_length = 25
        default_frag_length = 0 # replaced by read length

        self.c = cr.CountReadsPerBin([self.bamFile1, self.bamFile2],
                             binLength = bin_length,
                             defaultFragmentLength=default_frag_length,
                             stepSize=step_size)

    def test_countReadsInRegions_worker(self):
        step_size = 50
        bin_length = 25
        self.c.skipZeros = False
        resp = self.c.countReadsInRegions_worker(self.chrom, 0, 200)

        nt.assert_equal(resp, np.array([[ 0.,  0.],
                                        [ 0.,  1.],
                                        [ 1.,  1.],
                                        [ 1.,  2.]]))
        

        # step size = 200
        # bin length = 200
        # in other words, count the reads over the whole region
        # 2 for the first case, and 4 for the other
        self.c.stepSize = 200
        self.c.binLength = 200
        resp = self.c.countReadsInRegions_worker(self.chrom, 0, 200)
        nt.assert_equal(resp, np.array([[ 2., 4.]]))

    def test_countReadsInReginos_min_mapping_quality(self):
        # Test min mapping quality.
        self.c.minMappingQuality=40
        self.c.skipZeros = False

        resp = self.c.countReadsInRegions_worker(self. chrom, 0, 200)
        nt.assert_equal(resp, np.array([[ 0.,  0.,  0.,  1.],
                                        [ 0.,  0.,  0.,  1.]]).T)
    
    def test_countReadsInRegions_worker_ignore_duplicates(self):

        # Test ignore duplicates
        self.c.skipZeros = False
        self.c.ignoreDuplicates = True
        resp = self.c.countReadsInRegions_worker(self. chrom, 0, 200)

        nt.assert_equal(resp, np.array([[ 0.,  0.,  1.,  1.],
                                        [ 0.,  1.,  1.,  1.]]).T)
    
    def test_countReadsInRegions_worker_ignore_bed_regions(self):
        #Test bed regions:
        bed_regions = [(self. chrom, 10, 20), (self. chrom, 150, 160)]
        self.c.skipZeros = False
        resp =self.c.countReadsInRegions_worker(self. chrom, 0, 200, bed_regions_list=bed_regions)
        nt.assert_equal(resp, np.array([[ 0.,  1.],
                                        [ 0.,  2.]]).T)



