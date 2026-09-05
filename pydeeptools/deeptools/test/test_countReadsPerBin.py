# from unittest import TestCase

import deeptools.countReadsPerBin as cr
import numpy as np
import numpy.testing as nt
import os.path
import pytest

__author__ = 'Fidel'


@pytest.mark.parametrize("bc", ["bam", 'cram'])
class TestCountReadsPerBin():

    def ifiles(self, ext='bam'):
        root = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
        bamFile1 = root + "testA." + ext
        bamFile2 = root + "testB." + ext
        bamFile_PE = root + "test_paired2." + ext
        chrom = '3R'
        step_size = 50
        bin_length = 25
        c = cr.CountReadsPerBin(
            [bamFile1, bamFile2],
            binLength=bin_length,
            stepSize=step_size
        )
        return c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length
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

    def test_count_reads_in_region(self, bc):
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length = self.ifiles(bc)
        c.skipZeros = False
        resp, _ = c.count_reads_in_region(chrom, 0, 200)

        nt.assert_equal(resp, np.array([[0, 0.],
                                        [0, 1.],
                                        [1, 1.],
                                        [1, 2.]]))

    def test_count_reads_in_region_extension_1(self, bc):
        """
        In this case when read extension is smaller than read length
        extension is turned off and a warning is printed.
        """
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length = self.ifiles(bc)
        c = cr.CountReadsPerBin(
            [bamFile1, bamFile2],
            binLength=1,
            stepSize=50,
            extendReads=25
        )

        resp, _ = c.count_reads_in_region(chrom, 0, 200)

        nt.assert_equal(resp, np.array([[0, 0.],
                                        [0, 1.],
                                        [1, 1.],
                                        [1, 2.]]))

    def test_count_reads_in_region_total(self, bc):
        """ count the reads over the whole region
        2 for the first case, and 4 for the second
        """
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length = self.ifiles(bc)
        c.skipZeros = False
        c.stepSize = 200
        c.binLength = 200
        resp, _ = c.count_reads_in_region(chrom, 0, 200)
        nt.assert_equal(resp, np.array([[2, 4.]]))

    def test_countReadsInRegions_min_mapping_quality(self, bc):
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length = self.ifiles(bc)
        # Test min mapping quality.
        c.minMappingQuality = 40
        c.skipZeros = False

        resp, _ = c.count_reads_in_region(chrom, 0, 200)
        nt.assert_equal(resp, np.array([[0, 0, 0, 1.],
                                        [0, 0, 0, 1.]]).T)

    def test_count_reads_in_region_ignore_duplicates(self, bc):
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length = self.ifiles(bc)
        # Test ignore duplicates
        c.skipZeros = False
        c.ignoreDuplicates = True
        resp, _ = c.count_reads_in_region(chrom, 0, 200)

        nt.assert_equal(resp, np.array([[0, 0, 1, 1.],
                                        [0, 1, 1, 1.]]).T)

    def test_count_reads_in_region_ignore_bed_regions(self, bc):
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length = self.ifiles(bc)
        # Test bed regions:
        bed_regions = [[chrom, [(10, 20)], "."], [chrom, [(150, 160)], "."]]
        c.skipZeros = False
        c.binLength = 10
        resp, _ = c.count_reads_in_region(chrom, 0, 200, bed_regions_list=bed_regions)
        nt.assert_equal(resp, np.array([[0, 1.],
                                        [0, 2.]]).T)

    def test_get_coverage_of_region_sam_flag_include(self, bc):
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length = self.ifiles(bc)
        c.samFlag_include = 16  # include reverse reads only
        c.bamFilesList = [bamFile1]
        resp, _ = c.count_reads_in_region(chrom, 0, 200)
        nt.assert_array_equal(resp, np.array([[0], [0], [0], [1]]))

    def test_get_coverage_of_region_sam_flag_exclude(self, bc):
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length = self.ifiles(bc)
        c.samFlag_exclude = 16  # exclude reverse reads
        c.bamFilesList = [bamFile1]
        resp, _ = c.count_reads_in_region(chrom, 0, 200)
        nt.assert_array_equal(resp, np.array([[0], [0], [1], [0]]))

    def test_get_coverage_of_region_large_bin(self, bc):
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length = self.ifiles(bc)
        c.bamFilesList = [bamFile2]
        c.binLength = 200
        c.stepSize = 200
        resp, _ = c.count_reads_in_region(chrom, 0, 200)
        nt.assert_array_equal(resp, np.array([[4]]))

    def test_get_coverage_of_region_ignore_duplicates(self, bc):
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length = self.ifiles(bc)
        c.ignoreDuplicates = True
        c.bamFilesList = [bamFile2]
        resp, _ = c.count_reads_in_region(chrom, 0, 200)
        nt.assert_array_equal(resp, np.array([[0.],
                                              [1.],
                                              [1.],
                                              [1.]]))

        # check zero to nans
        c.zerosToNans = True
        resp, _ = c.count_reads_in_region(chrom, 0, 200)
        nt.assert_array_equal(resp, np.array([[np.nan],
                                              [1.],
                                              [1.],
                                              [1.]]))

    def test_get_coverage_of_region_split_read(self, bc):
        """
        The bamFile1 contains a read at position 10
        with the following CIGAR: 10S20M10N10M10S
        that maps to a chromosome named chr_cigar.
        """
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length = self.ifiles(bc)
        # turn of read extension
        c.extendPairedEnds = False
        c.bamFilesList = [bamFile1]
        c.binLength = 10
        c.stepSize = 10
        resp, _ = c.count_reads_in_region('chr_cigar', 0, 100)
        nt.assert_array_equal(resp, np.array([[0.],
                                              [1.],
                                              [1.],
                                              [0.],
                                              [1.],
                                              [0.],
                                              [0.],
                                              [0.],
                                              [0.],
                                              [0.]]))

    def test_get_coverage_of_region_zeros_to_nan(self, bc):
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length = self.ifiles(bc)
        c.zerosToNans = True
        resp, _ = c.count_reads_in_region(chrom, 0, 200)

        nt.assert_equal(resp, np.array([[np.nan, np.nan],
                                        [np.nan, 1],
                                        [1, 1],
                                        [1, 2]]))

    def test_bed_file(self, bc):
        c, bamFile1, bamFile2, bamFile_PE, chrom, step_size, bin_length = self.ifiles(bc)
        bed = "chr3R\t0\t10\nchr3R\t110\t120\nchr3R\t160\t180"
        import tempfile
        bed_file = tempfile.NamedTemporaryFile(suffix=".bed", delete=False, mode="w")
        bed_file.write(bed)
        bed_file.close()

        c = cr.CountReadsPerBin(
            [bamFile2],
            bedFile=[bed_file.name]
        )

        resp = c.run()
        nt.assert_equal(resp, np.array([[0.],
                                        [1.],
                                        [2.]]))

        import os
        os.unlink(bed_file.name)
