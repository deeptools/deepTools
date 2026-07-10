import pytest
from deeptools.computeMatrix2 import r_computematrix  # Adjust the import to your actual module
import os.path
import gzip
import hashlib

root = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
matrix = root + "computeMatrixOperations.mat.gz"
bed = root + "computeMatrixOperations.bed"
rbindMatrix1 = root + "somegenes.txt.gz"
rbindMatrix2 = root + "othergenes.txt.gz"
bigwig = root + "testA.bw"
outnpz = root + "output.mat.npz.gz"

def test_r_computematrix_referencePoint():
    mode = "reference-point"
    regionlis = [bed]
    bwlis = [bigwig]
    sampleslabel = ["sample1"]
    upstream = 1000
    downstream = 1000
    unscaled5prime = 0
    unscaled3prime = 0
    regionbodylength = 0
    binsize = 50
    missingdatazero = False
    metagene = False
    txnid = ""
    exonid = ""
    txniddesignator = ""
    scale = 1.0
    nanafterend = False
    skipzeros = False
    minthresh = 0.0
    maxthresh = 0.0
    averagetypebins = "mean"
    sortregions = "keep"
    sortusing = "mean"
    ortusingsamples = []
    referencepoint = "TSS"
    nproc = 1
    verbose = False
    ofile = outnpz

    result = r_computematrix(
        mode, regionlis, bwlis, sampleslabel, upstream, downstream, unscaled5prime, unscaled3prime,
        regionbodylength, binsize, missingdatazero, metagene, txnid, exonid, txniddesignator, scale,
        nanafterend, skipzeros, minthresh, maxthresh, averagetypebins, sortregions, sortusing,
        ortusingsamples, referencepoint, nproc, verbose, ofile
    )



    with gzip.open(ofile, 'rb') as f:
        file_content = f.read()
        h = hashlib.md5(file_content).hexdigest()

    expectedh = '4f1a2ce422d5b74fb6b75a81916929db'
    assert h == expectedh

    os.remove(ofile)

def test_r_computematrix_scale():
    mode = "scale-regions"
    regionlis = [bed]
    bwlis = [bigwig]
    sampleslabel = ["sample1"]
    upstream = 1000
    downstream = 1000
    unscaled5prime = 0
    unscaled3prime = 0
    regionbodylength = 0
    binsize = 50
    missingdatazero = False
    metagene = False
    txnid = ""
    exonid = ""
    txniddesignator = ""
    scale = 1.0
    nanafterend = False
    skipzeros = False
    minthresh = 0.0
    maxthresh = 0.0
    averagetypebins = "mean"
    sortregions = "keep"
    sortusing = "mean"
    ortusingsamples = []
    referencepoint = "TSS"
    nproc = 1
    verbose = False
    ofile = outnpz

    result = r_computematrix(
        mode, regionlis, bwlis, sampleslabel, upstream, downstream, unscaled5prime, unscaled3prime,
        regionbodylength, binsize, missingdatazero, metagene, txnid, exonid, txniddesignator, scale,
        nanafterend, skipzeros, minthresh, maxthresh, averagetypebins, sortregions, sortusing,
        ortusingsamples, referencepoint, nproc, verbose, ofile
    )

    with gzip.open(ofile, 'rb') as f:
        file_content = f.read()
        h = hashlib.md5(file_content).hexdigest()

    
    expectedh = '4f1a2ce422d5b74fb6b75a81916929db'
    assert h == expectedh

    os.remove(ofile)
    
