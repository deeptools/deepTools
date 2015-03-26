#from unittest import TestCase
from deeptools.countReadsPerBin import *
import numpy as np
import numpy.testing as nt

__author__ = 'Fidel'


class TestCountReadsInRegions_worker(object):
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

    def test_countReadsInRegions_worker(self):
        # step size = 50
        # bin length = 25
        resp = countReadsInRegions_worker(self.chrom, 0, 200,
                                          [self.bamFile1, self.bamFile2],
                                          50, 25, 0)

        nt.assert_equal(resp, np.array([[ 0.,  0.],
                                        [ 0.,  1.],
                                        [ 1.,  1.],
                                        [ 1.,  2.]]))
        

        # step size = 200
        # bin length = 200
        # in other words, count the reads over the whole region
        # 2 for the first case, and 4 for the other
        resp = countReadsInRegions_worker(self.chrom, 0, 200,
                                          [self. bamFile1, self. bamFile2], 200, 200, 0)
        nt.assert_equal(resp, np.array([[ 2., 4.]]))
    
        # Test min mapping quality.
        resp = countReadsInRegions_worker(self. chrom, 0, 200,
            [self. bamFile1, self. bamFile2], 50, 25, 0, minMappingQuality=40)
        nt.assert_equal(resp, np.array([[ 0.,  0.,  0.,  1.],
                                        [ 0.,  0.,  0.,  1.]]).T)
    
    def test_countReadsInRegions_worker_ignore_duplicates(self):

        # Test ignore duplicates
        # step size = 50
        # bin length = 25

        resp = countReadsInRegions_worker(self. chrom, 0, 200,
            [self. bamFile1, self. bamFile2], 50, 25, 0, ignoreDuplicates=True)

        nt.assert_equal(resp, np.array([[ 0.,  0.,  1.,  1.],
               [ 0.,  1.,  1.,  1.]]).T)
    
    def test_countReadsInRegions_worker_ignore_bed_regions(self):
        #Test bed regions:
        bedRegions = [(self. chrom, 10, 20), (self. chrom, 150, 160)]
        resp =countReadsInRegions_worker(self. chrom, 0, 200,
            [self. bamFile1, self. bamFile2], 0, 200, 0, bedRegions=bedRegions)
        nt.assert_equal(resp, np.array([[ 0.,  1.],
               [ 0.,  2.]]).T)



