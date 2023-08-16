from unittest import TestCase
import os
import pytest
import deeptools.writeBedGraph as wr
from deeptools.writeBedGraph import scaleCoverage

@pytest.mark.parametrize("bc", ["bam", 'cram'])
class TestWriteBedGraph():
    def ifiles(self, ext='bam'):
        root = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
        bamFile1 = root + "testA." + ext
        bamFile2 = root + "testB." + ext
        bamFile_PE = root + "test_paired2." + ext
        chrom = '3R'
        step_size = 50
        bin_length = 50
        func_args = {'scaleFactor': 1.0}
        c = wr.WriteBedGraph([bamFile1],
                                  binLength=bin_length,
                                  stepSize=step_size)
        return c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length, func_args

    def test_writeBedGraph_worker(self, bc):
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length, func_args = self.ifiles(bc)
        c.zerosToNans = False
        c.skipZeros = False

        tempFile = c.writeBedGraph_worker(chrom, 0, 200, scaleCoverage, func_args)
        _foo = open(tempFile[3], 'r')
        res = _foo.readlines()
        _foo.close()
        assert res == ['3R\t0\t100\t0\n', '3R\t100\t200\t1\n']
        os.remove(tempFile[3])

    def test_writeBedGraph_worker_zerotonan(self, bc):
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length, func_args = self.ifiles(bc)
        # turn on zeroToNan
        c.zerosToNans = True
        tempFile2 = c.writeBedGraph_worker(chrom, 0, 200, scaleCoverage, func_args)
        _foo = open(tempFile2[3], 'r')
        res = _foo.readlines()
        _foo.close()
        assert res == ['3R\t100\t200\t1\n']
        os.remove(tempFile2[3])

    def test_writeBedGraph_worker_scaling(self, bc):
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length, func_args = self.ifiles(bc)
        func_args = {'scaleFactor': 3.0}
        tempFile = c.writeBedGraph_worker(chrom, 0, 200, scaleCoverage, func_args)
        _foo = open(tempFile[3], 'r')
        res = _foo.readlines()
        _foo.close()
        assert res == ['3R\t0\t100\t0\n', '3R\t100\t200\t3\n']
        os.remove(tempFile[3])

    def test_writeBedGraph_worker_ignore_duplicates(self, bc):
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length, func_args = self.ifiles(bc)
        c = wr.WriteBedGraph([bamFile2],
                                  binLength=bin_length,
                                  stepSize=step_size, ignoreDuplicates=True)
        c.zerosToNans = True

        tempFile = c.writeBedGraph_worker(chrom, 0, 200, scaleCoverage, func_args)
        _foo = open(tempFile[3], 'r')
        res = _foo.readlines()
        _foo.close()
        assert res == ['3R\t50\t200\t1\n']
        os.remove(tempFile[3])

    def test_writeBedGraph_worker_smoothing(self, bc):
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length, func_args = self.ifiles(bc)
        c.binLength = 20
        c.stepSize = 20
        c.smoothLength = 60
        tempFile = c.writeBedGraph_worker(chrom, 100, 200, scaleCoverage, func_args)
        _foo = open(tempFile[3], 'r')
        res = _foo.readlines()
        _foo.close()
        assert res == ['3R\t100\t120\t1\n', '3R\t120\t180\t1.33333\n', '3R\t180\t200\t1\n']
        os.remove(tempFile[3])

    def test_writeBedGraph_cigar(self, bc):
        """
        The bamFile1 contains a read at position 10
        with the following CIGAR: 10S20M10N10M10S
        that maps to a chromosome named chr_cigar.
        """
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length, func_args = self.ifiles(bc)
        # turn of read extension
        c.extendPairedEnds = False
        c.binLength = 10
        c.stepSize = 10
        tempFile = c.writeBedGraph_worker('chr_cigar', 0, 100, scaleCoverage, func_args)
        _foo = open(tempFile[3], 'r')
        res = _foo.readlines()
        _foo.close()

        # the sigle read is split into bin 10-30, and then 40-50
        assert res == ['chr_cigar\t0\t10\t0\n',
                           'chr_cigar\t10\t30\t1\n',
                           'chr_cigar\t30\t40\t0\n',
                           'chr_cigar\t40\t50\t1\n',
                           'chr_cigar\t50\t100\t0\n']
        os.remove(tempFile[3])
