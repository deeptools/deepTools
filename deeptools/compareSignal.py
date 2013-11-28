#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
import multiprocessing

# own tools
from deeptools import writeBedGraph, bamHandler
from deeptools.countReadsPerBin import getNumReadsPerBin
from scipy.stats import poisson

debug = 1

"""
 NUM = ChIP treatment 
 DEN = Input

 calculate p-values assuming 
      NUM = signal-N + noise-N (where contribution 
                                   of signal-N & noise-N is unknown)

      DEN = signal-D + noise-D (where contribution of signal-D & noise-D is unknown)

      p-value = probability of null hypothesis = signal-N is the same as signal-D

          1. identify regions of NUM with highest values -> greatest contribution of signal-N (so least contribution of noise-N)
          2. if signal-N = signal-D, these should also correspond to regions of DEN with greatest contribution of signal-D
               (even if contribution of signal-D to DEN is still small compared to noise-D)
          3. calculate mean(NUM) at these regions = mean(noise-N) + mean(high signal-N) 
                                                  = mean(noise-N) + mean(signal-N) + mean(extra signal-N at these regions)
                                                  = mean(NUM) + mean(extra signal-N at these regions)
          4. calculate mean(DEN) at these regions = mean(noise-D) + mean(high signal-D) 
                                                  = mean(noise-D) + mean(signal-D) + mean(extra signal-D at these regions)
                                                  = mean(DEN) + mean(extra signal-D at these regions)
          5. estimate ratio of magnitudes of signal in NUM & DEN:
               since magnitude of extra signal at particular regions is proportional to magnitude of mean signal,
               signalRatio = mean(signal-N)/mean(signal-D)
                           = mean(extra signal-N at these regions)/mean(extra signal-D at these regions)
          6. now, for every genomic interval, calculate lambda = expected DEN signal, from NUM signal based on null hypothesis:
               lambda = expected DEN signal = mean(DEN) + (NUM signal - mean(NUM))/signalRatio
          7. use calculated lambda to determine p-value that signal ≤ observed DEN signal

SYNOPSIS: compareSignal [options] DENOMINATORFILE.bed NUMERATORFILE.bed

OPTIONS:   -c = input files sorted by CHROM (requires less memory)
           -f = FRAGMENT size = 2x offset for shifting tag positions (default = 300)
                (set to 'pe' for paired-end data)
           -s = STEP size to write output (default =50)
           -w = WINDOW size to count tags (default = 200), or
           -r = REGION file of genomic intervals to count tags
           -m = MAXIMUM ratio of underlying NUM/DEN signal (default = 10)
use 'stdin' to read either file from STDIN

NOTE: retains duplicate BED values (=same chrom, chromStart values)
"""

def computeLambda(tileCoverage, args):
    """
    This function is called by the writeBedGraph workers for every 
    tile in the genome that is considered
    """

    treatmentWindowTags = tileCoverage[0]
    controlWindowTags = tileCoverage[1]
    treatmentExtraSignalTags = treatmentWindowTags - args['treatmentMean']
    
    controlLambda = args['controlMean'] + (treatmentExtraSignalTags * args['controlSignalRatio'])
    
    return controlLambda


def computePvalue(tileCoverage, args):
    """
    This function is called by the writeBedGraph workers for every 
    tile in the genome that is considered
    """
#    if tileCoverage == (0,0):
#        return np.nan

    treatmentWindowTags = tileCoverage[0]
    controlWindowTags = tileCoverage[1]
    if controlWindowTags == 0:
        return np.nan

    treatmentExtraSignalTags = treatmentWindowTags - args['treatmentMean']
    
    controlLambda = args['controlMean'] + (treatmentExtraSignalTags * args['controlSignalRatio'])
    
    log10pvalue = -1* poisson.logcdf(controlWindowTags, controlLambda) / np.log(10)
#    log10pvalue = -1* poisson.logsf(controlWindowTags, controlLambda) / np.log(10)

    return log10pvalue


def compareSignal(bamFilesList, binLength, numberOfSamples, defaultFragmentLength, 
                  outFileName, outFileFormat, outFileNameLambda=None, region=None,
                  extendPairedEnds=True,
                  numberOfProcessors=1, Nsigmas = 2, maxSignalRatio=10, verbose=False):
    
    bam1 = bamHandler.openBam(bamFilesList[0])
    genomeSize = sum(bam1.lengths)

    bam2 = bamHandler.openBam(bamFilesList[1])

    treatmentMapped = bam1.mapped
    controlMapped  =  bam2.mapped
    treatmentControlRatioMapped = float(treatmentMapped) / controlMapped

    # 1. Get a table containing number of reads in a sample from the genome.
    #    Only regions for which both samples have non zero counts are considered

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

    print "Treatment mean: {:.2f}, Treatment total:{:.2f}".format(treatmentMean, treatmentTotal)
    print "Control mean: {:.2f}, Control total:{}".format(controlMean, controlTotal)
    print "the ratio of treatment vs. control for enriched regions is: {:.2f}".format(treatmentSignalRatio)
    print "the ratio of treatment vs. control ratio: {:.2f} (if based on mapped reads: {:.2f})".format(treatmentControlRatio, treatmentControlRatioMapped)

    

    funcArgs = {'controlMean': controlMean,
                'treatmentMean': treatmentMean,
                'controlSignalRatio': controlSignalRatio,
                'controlRatio': controlRatio,
                'treatmentControlRatio': treatmentControlRatio
                }


    writeBedGraph.writeBedGraph( bamFilesList,
                                 outFileName,
                                 defaultFragmentLength, computePvalue, 
                                 funcArgs, tileSize=binLength, region=region,
                                 format=outFileFormat,
                                 zerosToNans = False,
                                 numberOfProcessors=numberOfProcessors,
                                 extendPairedEnds=extendPairedEnds)

    if outFileNameLambda:
        writeBedGraph.writeBedGraph( bamFilesList,
                                     outFileNameLambda,
                                     defaultFragmentLength, computeLambda, 
                                     funcArgs, tileSize=binLength, region=region,
                                     format=outFileFormat,
                                     zerosToNans = False,
                                     numberOfProcessors=numberOfProcessors,
                                     extendPairedEnds=extendPairedEnds)


