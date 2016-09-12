#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np

# own packages
from deeptools import bamHandler
import deeptools.countReadsPerBin as countR

old_settings = np.seterr(all='ignore')
debug = 0


def estimateScaleFactor(bamFilesList, binLength, numberOfSamples,
                        normalizationLength,
                        avg_method='median', blackListFileName=None, numberOfProcessors=1,
                        verbose=False, chrsToSkip=[]):
    r"""
    Subdivides the genome into chunks to be analyzed in parallel
    using several processors. The code handles the creation of
    workers that compute fragment counts (coverage) for different
    regions and then collect and integrates the results.

    Parameters
    ----------
    bamFilesList : list
        list of bam files to normalize
    binLength : int
        the window size in bp, where reads are going to be
        counted.
    numberOfSamples : int
        number of sites to sample from the genome. For more info see
        the documentation of the CountReadsPerBin class
    normalizationLength : int
        length, in bp, to normalize the data.
        For a value of 1, on average
        1 read per base pair is found
    avg_method : str
        defines how the different values are to be summarized.
        The options are 'mean' and 'median'
    chrsToSkip : list
        name of the chromosomes to be excluded from the
        scale estimation. Usually the chrX is included.
    blackListFileName : str
        BED file containing blacklisted regions

    Returns
    -------
    dict
        Dictionary with the following keys::
            'size_factors'
            'size_factors_based_on_mapped_reads'
            'size_factors_SES'
            'size_factors_based_on_mean'
            'size_factors_based_on_median'
            'mean'
            'meanSES'
            'median'
            'reads_per_bin'
            'std'
            'sites_sampled'


    Examples
    --------
    >>> test = Tester()
    >>> bin_length = 50
    >>> num_samples = 4
    >>> _dict = estimateScaleFactor([test.bamFile1, test.bamFile2], bin_length, num_samples,  1)
    >>> _dict['size_factors']
    array([ 1. ,  0.5])
    >>> _dict['size_factors_based_on_mean']
    array([ 1. ,  0.5])
    """

    assert len(bamFilesList) == 2, "SES scale factors are only defined for 2 files"

    bamFilesHandlers = [bamHandler.openBam(x) for x in bamFilesList]
    mappedReads = [x.mapped for x in bamFilesHandlers]

    sizeFactorBasedOnMappedReads = np.array(mappedReads, dtype='float64')

    sizeFactorBasedOnMappedReads = sizeFactorBasedOnMappedReads.min() / sizeFactorBasedOnMappedReads

    cr = countR.CountReadsPerBin(bamFilesList,
                                 binLength=binLength,
                                 numberOfSamples=numberOfSamples,
                                 extendReads=False,
                                 blackListFileName=blackListFileName,
                                 numberOfProcessors=numberOfProcessors,
                                 verbose=verbose,
                                 chrsToSkip=chrsToSkip)

    try:
        num_reads_per_bin = cr.run()
    except Exception as detail:
        exit("*ERROR*: {}".format(detail))

    sitesSampled = len(num_reads_per_bin)

    # the transpose is taken to easily iterate by columns which are now
    # converted to rows
    num_reads_per_bin = num_reads_per_bin.transpose()
#    np.savetxt("/home/ramirez/tmp/test.num_reads", num_reads_per_bin)
    # size factors based on order statistics
    # see Signal extraction scaling (SES) method in: Diaz et al (2012)
    # Normalization, bias correction, and peak calling for ChIP-seq.
    # Statistical applications in genetics and molecular biology, 11(3).

    # using the same names as in Diaz paper
    # p refers to ChIP, q to input

    p = np.sort(num_reads_per_bin[0, :]).cumsum()
    q = np.sort(num_reads_per_bin[1, :]).cumsum()

    # p[-1] and q[-1] are the maximum values in the  arrays.
    # both p and q are normalized by this value
    diff = np.abs(p / p[-1] - q / q[-1])
    # get the lowest rank for wich the difference is the maximum
    maxIndex = np.flatnonzero(diff == diff.max())[0]
    # Take a lower rank to move to a region with probably
    # less peaks and more background.
    maxIndex = int(maxIndex * 0.8)
    while(maxIndex < len(p)):
        # in rare cases the maxIndex maps to a zero value.
        # In such cases, the next index is used until
        # a non zero value appears.
        cumSum = np.array([float(p[maxIndex]), float(q[maxIndex])])
        if cumSum.min() > 0:
            break
        maxIndex += 1

    meanSES = [np.mean(np.sort(num_reads_per_bin[0, :])[:maxIndex]),
               np.mean(np.sort(num_reads_per_bin[1, :])[:maxIndex])]

    # the maxIndex may be too close to the the signal regions
    # so i take a more conservative approach by taking a close number

    sizeFactorsSES = cumSum.min() / cumSum
    median = np.median(num_reads_per_bin, axis=1)

    # consider only those read numbers that are below the 90
    # percentile to stimate the
    # mean and std
    mean = []
    std = []
    for values in num_reads_per_bin:
        maxNumReads = (np.percentile(values, 90))
        if maxNumReads == 0:
            maxNumReads = (np.percentile(values, 99))
            if maxNumReads == 0:
                print("all genomic regions sampled from one ")
                "of the bam files have no reads.\n"
                values = values[values <= maxNumReads]

        mean.append(np.mean(values))
        std.append(np.std(values))

    mean = np.array(mean)
    readsPerBin = mean if avg_method == 'mean' else median

    if min(median) == 0:
        idx_zero = [ix + 1 for ix, value in enumerate(median) if value == 0]
        exit("\n*ERROR*: The median coverage computed is zero for sample(s) #{}\n"
             "Try selecting a larger sample size or a region with coverage\n".format(idx_zero))

    sizeFactor = sizeFactorsSES
    return {'size_factors': sizeFactor,
            'size_factors_based_on_mapped_reads': sizeFactorBasedOnMappedReads,
            'size_factors_SES': sizeFactorsSES,
            'size_factors_based_on_mean': mean.min() / mean,
            'size_factors_based_on_median': median.min() / median,
            'mean': mean,
            'meanSES': meanSES,
            'median': median,
            'reads_per_bin': readsPerBin,
            'std': std,
            'sites_sampled': sitesSampled}


class Tester(object):

    def __init__(self):
        self.root = os.path.dirname(os.path.abspath(__file__)) + "/test/test_data/"
        self.bamFile1 = self.root + "testA.bam"
        self.bamFile2 = self.root + "testB.bam"
        global debug
        debug = 0
        self.chrom = '3R'
