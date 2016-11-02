# from unittest import TestCase

import deeptools.computeMatrixOperations as cmo
import os
import hashlib
import gzip
import json

__author__ = 'Devon'

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"


def getHeader(fp):
    s = fp.readline()
    if isinstance(s, bytes):
        s = s.decode()
    s = s[1:]
    return json.loads(s)


class TestComputeMatrixOperations(object):
    def setUp(self):
        self.root = ROOT
        self.matrix = self.root + "computeMatrixOperations.mat.gz"
        self.bed = self.root + "computeMatrixOperations.bed"

    def testSubset(self):
        """
        computeMatrixOperations subset
        """
        dCorrect = {"verbose": True, "scale": 1, "skip zeros": False, "nan after end": False, "sort using": "mean", "unscaled 5 prime": 0, "body": 1000, "sample_labels": ["SRR648667.forward", "SRR648668.forward", "SRR648669.forward", "SRR648670.forward"], "downstream": 0, "unscaled 3 prime": 0, "group_labels": ["genes"], "bin size": 10, "upstream": 0, "group_boundaries": [0, 196], "sample_boundaries": [0, 100, 200, 300, 400], "max threshold": None, "ref point": None, "min threshold": None, "sort regions": "no", "proc number": 20, "bin avg type": "mean", "missing data as zero": False}
        oname = "/tmp/subset.mat.gz"
        args = "subset -m {} --sample SRR648667.forward SRR648668.forward SRR648669.forward SRR648670.forward -o {}".format(self.matrix, oname)
        args = args.split()
        cmo.main(args)
        f = gzip.GzipFile(oname)
        d = getHeader(f)  # Skip the header, which can be in a different order
        h = hashlib.md5(f.read()).hexdigest()
        f.close()
        assert(d == dCorrect)
        assert(h == "edb3c8506c3f27ebb8c7ddf94d5ba594")
        os.remove(oname)

    def testfilterStrand(self):
        """
        computeMatrixOperations filterStrand
        """
        dCorrect = {"verbose": True, "scale": 1, "skip zeros": False, "nan after end": False, "sort using": "mean", "unscaled 5 prime": 0, "body": 1000, "sample_labels": ["SRR648667.forward", "SRR648668.forward", "SRR648669.forward", "SRR648670.forward", "SRR648667.reverse", "SRR648668.reverse", "SRR648669.reverse", "SRR648670.reverse"], "downstream": 0, "unscaled 3 prime": 0, "group_labels": ["genes"], "bin size": 10, "upstream": 0, "group_boundaries": [0, 107], "sample_boundaries": [0, 100, 200, 300, 400, 500, 600, 700, 800], "max threshold": None, "ref point": None, "min threshold": None, "sort regions": "no", "proc number": 20, "bin avg type": "mean", "missing data as zero": False}
        oname = "/tmp/filterStrand1.mat.gz"
        args = "filterStrand -m {} -o {} --strand +".format(self.matrix, oname)
        args = args.split()
        cmo.main(args)
        f = gzip.GzipFile(oname)
        d = getHeader(f)  # Skip the header, which can be in a different order
        h = hashlib.md5(f.read()).hexdigest()
        f.close()
        assert(d == dCorrect)
        assert(h == "300f8000be5b5f51e803b57ef08f1c9e")
        os.remove(oname)

        dCorrect = {u'verbose': True, u'scale': 1, u'skip zeros': False, u'nan after end': False, u'sort using': u'mean', u'unscaled 5 prime': 0, u'body': 1000, u'sample_labels': [u'SRR648667.forward', u'SRR648668.forward', u'SRR648669.forward', u'SRR648670.forward', u'SRR648667.reverse', u'SRR648668.reverse', u'SRR648669.reverse', u'SRR648670.reverse'], u'downstream': 0, u'unscaled 3 prime': 0, u'group_labels': [u'genes'], u'bin size': 10, u'upstream': 0, u'group_boundaries': [0, 89], u'sample_boundaries': [0, 100, 200, 300, 400, 500, 600, 700, 800], u'missing data as zero': False, u'ref point': None, u'min threshold': None, u'sort regions': u'no', u'proc number': 20, u'bin avg type': u'mean', u'max threshold': None}
        oname = "/tmp/filterStrand2.mat.gz"
        args = "filterStrand -m {} -o {} --strand -".format(self.matrix, oname)
        args = args.split()
        cmo.main(args)
        f = gzip.GzipFile(oname)
        d = getHeader(f)  # Skip the header, which can be in a different order
        h = hashlib.md5(f.read()).hexdigest()
        f.close()
        assert(d == dCorrect)
        assert(h == "0a6ca070a5ba4564f1ab950ac3b7c8f1")
        os.remove(oname)

    def testrbind(self):
        """
        computeMatrixOperations rbind
        """
        dCorrect = {"verbose": True, "scale": 1, "skip zeros": False, "nan after end": False, "sort using": "mean", "unscaled 5 prime": 0, "body": 1000, "sample_labels": ["SRR648667.forward", "SRR648668.forward", "SRR648669.forward", "SRR648670.forward", "SRR648667.reverse", "SRR648668.reverse", "SRR648669.reverse", "SRR648670.reverse"], "downstream": 0, "unscaled 3 prime": 0, "group_labels": ["genes"], "bin size": 10, "upstream": 0, "group_boundaries": [0, 392], "sample_boundaries": [0, 100, 200, 300, 400, 500, 600, 700, 800], "max threshold": None, "ref point": None, "min threshold": None, "sort regions": "no", "proc number": 20, "bin avg type": "mean", "missing data as zero": False}
        oname = "/tmp/rbind.mat.gz"
        args = "rbind -m {0} {0} -o {1}".format(self.matrix, oname)
        args = args.split()
        cmo.main(args)
        f = gzip.GzipFile(oname)
        d = getHeader(f)  # Skip the header, which can be in a different order
        h = hashlib.md5(f.read()).hexdigest()
        f.close()
        assert(d == dCorrect)
        assert(h == "3dd96c7b05e0ca5ada21212defe57fba")
        os.remove(oname)

    def testcbind(self):
        """
        computeMatrixOperations cbind
        """
        dCorrect = {"verbose": True, "scale": 1, "skip zeros": False, "nan after end": False, "sort using": "mean", "unscaled 5 prime": 0, "body": 1000, "sample_labels": ["SRR648667.forward", "SRR648668.forward", "SRR648669.forward", "SRR648670.forward", "SRR648667.reverse", "SRR648668.reverse", "SRR648669.reverse", "SRR648670.reverse", "SRR648667.forward", "SRR648668.forward", "SRR648669.forward", "SRR648670.forward", "SRR648667.reverse", "SRR648668.reverse", "SRR648669.reverse", "SRR648670.reverse"], "downstream": 0, "unscaled 3 prime": 0, "group_labels": ["genes"], "bin size": 10, "upstream": 0, "group_boundaries": [0, 196], "sample_boundaries": [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600], "max threshold": None, "ref point": None, "min threshold": None, "sort regions": "no", "proc number": 20, "bin avg type": "mean", "missing data as zero": False}
        oname = "/tmp/filterStrand.mat.gz"
        args = "cbind -m {0} {0} -o {1}".format(self.matrix, oname)
        args = args.split()
        cmo.main(args)
        f = gzip.GzipFile(oname)
        d = getHeader(f)  # Skip the header, which can be in a different order
        h = hashlib.md5(f.read()).hexdigest()
        f.close()
        assert(d == dCorrect)
        assert(h == "e55d89704bb16a11f366663a8fd90a47")
        os.remove(oname)

    def testsort(self):
        """
        computeMatrixOperations sort
        """
        dCorrect = {"verbose": True, "scale": 1, "skip zeros": False, "nan after end": False, "sort using": "mean", "unscaled 5 prime": 0, "body": 1000, "sample_labels": ["SRR648667.forward", "SRR648668.forward", "SRR648669.forward", "SRR648670.forward", "SRR648667.reverse", "SRR648668.reverse", "SRR648669.reverse", "SRR648670.reverse"], "downstream": 0, "unscaled 3 prime": 0, "group_labels": ["genes"], "bin size": 10, "upstream": 0, "group_boundaries": [0, 196], "sample_boundaries": [0, 100, 200, 300, 400, 500, 600, 700, 800], "max threshold": None, "ref point": None, "min threshold": None, "sort regions": "no", "proc number": 20, "bin avg type": "mean", "missing data as zero": False}
        oname = "/tmp/sorted.mat.gz"
        args = "sort -m {} -o {} -R {}".format(self.matrix, oname, self.bed)
        args = args.split()
        cmo.main(args)
        f = gzip.GzipFile(oname)
        d = getHeader(f)  # Skip the header, which can be in a different order
        h = hashlib.md5(f.read()).hexdigest()
        f.close()
        assert(d == dCorrect)
        assert(h == "10ea07d1aa58f44625abe2142ef76094")
        os.remove(oname)
