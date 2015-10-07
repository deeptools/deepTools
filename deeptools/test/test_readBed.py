import deeptools.readBed as readbed
from types import *


class TestReadBed(object):
    def setUp(self):
        self.root = "./test/test_data/"
        self.bed3_file  = self.root + "test.bed3"

        self.bed = readbed.ReadBed(open(self.bed3_file))
    def test_bed3(self):
        """
        in case of bed3 file, the columns 4-6 are added
        :return:
        """
        interval = self.bed.next()
        #import ipdb;ipdb.set_trace()
        # check the original file line
        assert interval.line == 'chr1\t1\t10\n'

        # check that the 4-6 columns are added
        assert interval.chrom == "chr1"
        assert interval.start ==  1
        assert interval.end == 10
        assert interval.name ==  "."
        assert interval.score == 0
        assert interval.strand == "."

    def test_casting(self):
        # chromosome name should always be a strand
        fields = self.bed.get_values_from_line("1\t0\t10\tNAME\t0.0\t+\n")
        assert type(fields['chrom']) == StringType
        #import ipdb;ipdb.set_trace()
        assert type(fields['start']) == IntType
        assert type(fields['end']) == IntType
        assert type(fields['name']) == StringType
        assert type(fields['score']) == FloatType
        assert type(fields['strand']) == StringType

        # in this test, the name and the strand should be cast
        # as strings
        fields = self.bed.get_values_from_line("1\t0\t10\t10\t0.0\t0\n")
        assert type(fields['name']) == StringType
        assert type(fields['strand']) == StringType

        # test skip of invalid field for start position
        fields = self.bed.get_values_from_line("1\twrong\t10\t10\t0.0\t0\n")
        assert fields == dict()
