from unittest import TestCase
from nose.tools import *
import os

import deeptools.writeBedGraph as wr
from deeptools.writeBedGraph import scaleCoverage

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"

__author__ = 'fidel'


class TestWriteBedGraph(TestCase):

    def setUp(self):
        """
        The distribution of reads between the two bam files is as follows.

        They cover 200 bp::

              0                              100                           200
              |------------------------------------------------------------|
            A                                ==============>
                                                            <==============


            B                 <==============               ==============>
                                             ==============>
                                                            ==============>
        """

        self.root = ROOT
        self.bamFile1 = self.root + "testA.bam"
        self.bamFile2 = self.root + "testB.bam"
        self.bamFile_PE = self.root + "test_paired2.bam"
        self.chrom = '3R'

        self.step_size = 50
        self.bin_length = 50
        self.func_args = {'scaleFactor': 1.0}

        self.c = wr.WriteBedGraph([self.bamFile1],
                                  binLength=self.bin_length,
                                  stepSize=self.step_size)

    def test_writeBedGraph_worker(self):
        self.c.zerosToNans = False
        self.c.skipZeros = False

        tempFile = self.c.writeBedGraph_worker('3R', 0, 200, scaleCoverage, self.func_args)
        _foo = open(tempFile, 'r')
        res = _foo.readlines()
        _foo.close()
        assert_equal(res, ['3R\t0\t100\t0.00\n', '3R\t100\t200\t1.00\n'])
        os.remove(tempFile)

    def test_writeBedGraph_worker_zerotonan(self):
        # turn on zeroToNan
        self.c.zerosToNans = True
        tempFile2 = self.c.writeBedGraph_worker('3R', 0, 200, scaleCoverage, self.func_args)
        _foo = open(tempFile2, 'r')
        res = _foo.readlines()
        _foo.close()
        assert_equal(res, ['3R\t100\t200\t1.00\n'])
        os.remove(tempFile2)

    def test_writeBedGraph_worker_scaling(self):
        func_args = {'scaleFactor': 3.0}
        tempFile = self.c.writeBedGraph_worker('3R', 0, 200, scaleCoverage, func_args)
        _foo = open(tempFile, 'r')
        res = _foo.readlines()
        _foo.close()
        assert_equal(res, ['3R\t0\t100\t0.00\n', '3R\t100\t200\t3.00\n'])
        os.remove(tempFile)

    def test_writeBedGraph_worker_ignore_duplicates(self):
        self.c = wr.WriteBedGraph([self.bamFile2],
                                  binLength=self.bin_length,
                                  stepSize=self.step_size, ignoreDuplicates=True)
        self.c.zerosToNans = True

        tempFile = self.c.writeBedGraph_worker('3R', 0, 200, scaleCoverage, self.func_args)
        _foo = open(tempFile, 'r')
        res = _foo.readlines()
        _foo.close()
        assert_equal(res, ['3R\t50\t200\t1.00\n'])
        os.remove(tempFile)

    def test_writeBedGraph_worker_smoothing(self):
        self.c.binLength = 20
        self.c.stepSize = 20
        self.c.smoothLength = 60
        tempFile = self.c.writeBedGraph_worker('3R', 100, 200, scaleCoverage, self.func_args)
        _foo = open(tempFile, 'r')
        res = _foo.readlines()
        _foo.close()
        assert_equal(res, ['3R\t100\t120\t1.00\n', '3R\t120\t180\t1.33\n', '3R\t180\t200\t1.00\n'])
        os.remove(tempFile)

    def test_writeBedGraph_cigar(self):
        """
        The bamFile1 contains a read at position 10
        with the following CIGAR: 10S20M10N10M10S
        that maps to a chromosome named chr_cigar.
        """

        # turn of read extension
        self.c.extendPairedEnds = False
        self.c.binLength = 10
        self.c.stepSize = 10
        tempFile = self.c.writeBedGraph_worker('chr_cigar', 0, 100, scaleCoverage, self.func_args)
        _foo = open(tempFile, 'r')
        res = _foo.readlines()
        _foo.close()

        # the sigle read is split into bin 10-30, and then 40-50
        assert_equal(res, ['chr_cigar\t0\t10\t0.00\n',
                           'chr_cigar\t10\t30\t1.00\n',
                           'chr_cigar\t30\t40\t0.00\n',
                           'chr_cigar\t40\t50\t1.00\n',
                           'chr_cigar\t50\t100\t0.00\n'])
        os.remove(tempFile)
