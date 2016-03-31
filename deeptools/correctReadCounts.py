#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import poisson
from deeptools import writeBedGraph
from deeptools.SES_scaleFactor import estimateScaleFactor

debug = 1
old_settings = np.seterr(all='ignore')


def computeLambda(tileCoverage, args):
    """
    This function is called by the writeBedGraph workers for every
    tile in the genome that is considered
    """

    treatmentWindowTags = tileCoverage[0]
    controlWindowTags = tileCoverage[1]
    treatmentExtraSignalTags = treatmentWindowTags - args['treatmentMean']

    controlLambda = args['controlMean'] + (treatmentExtraSignalTags * args['controlSignalRatio'])

    log10pvalue = -1 * poisson.logsf(controlWindowTags, controlLambda) / np.log(10)

    return log10pvalue


def computePvalue(tileCoverage, args):
    """
    This function is called by the writeBedGraph workers for every
    tile in the genome that is considered

    It computes a pvalue based on an expected lambda comming from
    the correction of treatment when the input is considered.
    """

    treatmentWindowTags = tileCoverage[0]
    controlWindowTags = tileCoverage[1]

    treatmentLambda = controlWindowTags * args['treatmentControlRatio']

    log10pvalue = min(300, -1 * poisson.logsf(treatmentWindowTags, treatmentLambda) / np.log(10))

    return log10pvalue


def computeCorrectedReadcounts(tileCoverage, args):
    """
    This function is called by the writeBedGraph workers for every
    tile in the genome that is considered

    It computes a pvalue based on an expected lambda comming from
    the correction of treatment when the input is considered.
    """

    treatmentWindowTags = tileCoverage[0]
    controlWindowTags = tileCoverage[1]
    treatmentCorrectedTags = treatmentWindowTags - args['treatmentControlRatio'] * (controlWindowTags - args['controlMean'])

    return treatmentCorrectedTags


def correctReadCounts(bamFilesList, binLength, numberOfSamples, defaultFragmentLength,
                      outFileName, outFileFormat, outFileNameCorr=None, region=None,
                      extendPairedEnds=True,
                      numberOfProcessors=1, Nsigmas=2, maxSignalRatio=10, blackListFileName=None, verbose=False):

    bam1 = writeBedGraph.openBam(bamFilesList[0])
    bam2 = writeBedGraph.openBam(bamFilesList[1])

    treatmentMapped = bam1.mapped
    controlMapped = bam2.mapped
    treatmentControlRatioMapped = float(treatmentMapped) / controlMapped

    # 1. Get a table containing number of reads in a sample from the genome.
    #    Only regions for which both samples have more than zero counts are considered

    scaleFactorsDict = estimateScaleFactor(bamFilesList,
                                           binLength,
                                           numberOfSamples,
                                           defaultFragmentLength,
                                           1,
                                           blackListFileName=blackListFileName,
                                           numberOfProcessors=numberOfProcessors,
                                           verbose=verbose)

    """
    num_reads_per_region = getNumReadsPerBin(bamFilesList, binLength, numberOfSamples, defaultFragmentLength, numberOfProcessors, skipZeros=True, verbose=verbose)
    if verbose:
        print "number of non-zero regions sampled: {}".format(num_reads_per_region.shape[0])

    # 2. get Mean and std of treatment (col1) and control (col2)

    treatmentMean, controlMean = np.mean(num_reads_per_region, axis=0) # axis=0: that measn by column
    treatmentStd, controlStd   = np.std(num_reads_per_region, axis=0)
    treatmentTotal, controlTotal   = np.sum(num_reads_per_region, axis=0)

    # 3. Calculate residual in treatment & control data, at regions for which treatment
    #    signal exceeds mean + std * Nsigmas
    #    (these are expected to be the regions at which the signal > mean-signal,
    #    so the residual signal is positive)

    overRows = np.where(num_reads_per_region[:,0].copy() >= treatmentMean + treatmentStd*Nsigmas )[0]
    over_Nsigma_regions = num_reads_per_region[overRows, :]

    treatmentSigMean, controlSigMean = np.mean(over_Nsigma_regions, axis=0)

    treatmentExtraSignal = treatmentSigMean - treatmentMean
    controlExtraSignal   = controlSigMean - controlMean

    treatmentControlRatio = float(treatmentTotal) / controlTotal
    adjSignalRatio = maxSignalRatio * treatmentControlRatio;
    treatmentSignalRatio = float(treatmentExtraSignal) / controlExtraSignal

    if treatmentSignalRatio < adjSignalRatio and treatmentSignalRatio > 0:
        treatmentSignalRatio = adjSignalRatio

    if treatmentSignalRatio < 1:
        raise NameError("estimated signal in control file {} is greater than estimated signal in treatmant file {}. Perhaps the file names are swapped?".format(bamFilesList[0], bamFilesList[1]))

    else:
        controlSignalRatio = 1.0/treatmentSignalRatio

    controlRatio = 1.0 / treatmentControlRatio

    """

#    scaleFactors = scaleFactorsDict['size_factors']

    treatmentMean, controlMean = scaleFactorsDict['meanSES']
    treatmentControlRatio = scaleFactorsDict['size_factors'][1] / scaleFactorsDict['size_factors'][0]
    treatmentSignalRatio = treatmentControlRatio
    controlRatio = controlSignalRatio = 1.0 / treatmentControlRatio
    treatmentTotal = treatmentMapped
    controlTotal = controlMapped

    print("Treatment mean: {:.2f}, Treatment total:{:.2f}".format(treatmentMean, treatmentTotal))
    print("Control mean: {:.2f}, Control total:{}".format(controlMean, controlTotal))
    print("the ratio of treatment vs. control for enriched regions is: {:.2f}".format(treatmentSignalRatio))
    print("the ratio of treatment vs. control ratio: {:.2f} (if based on mapped reads: {:.2f})".format(treatmentControlRatio, treatmentControlRatioMapped))

    funcArgs = {'controlMean': controlMean,
                'treatmentMean': treatmentMean,
                'controlSignalRatio': controlSignalRatio,
                'controlRatio': controlRatio,
                'treatmentControlRatio': treatmentControlRatio
                }

    writeBedGraph.writeBedGraph(bamFilesList,
                                outFileName,
                                defaultFragmentLength, computePvalue,
                                funcArgs, tileSize=binLength, region=region,
                                format=outFileFormat,
                                zerosToNans=False,
                                blackListFileName=blackListFileName,
                                numberOfProcessors=numberOfProcessors,
                                extendPairedEnds=extendPairedEnds)

    if outFileNameCorr:
        writeBedGraph.writeBedGraph(bamFilesList,
                                    outFileNameCorr,
                                    defaultFragmentLength, computeCorrectedReadcounts,
                                    funcArgs, tileSize=binLength, region=region,
                                    format=outFileFormat,
                                    zerosToNans=False,
                                    blackListFileName=blackListFileName,
                                    numberOfProcessors=numberOfProcessors,
                                    extendPairedEnds=extendPairedEnds)


def controlLambda(tileCoverage, args):

    controlLambda = args['controlMean'] + (tileCoverage[0] * args['controlSignalRatio'])
    return controlLambda
