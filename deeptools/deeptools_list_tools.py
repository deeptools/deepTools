#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
from deeptools._version import __version__


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
deepTools is a suite of python tools particularly developed for the efficient analysis of
high-throughput sequencing data, such as ChIP-seq, RNA-seq or MNase-seq.

Each tool should be called by its own name as in the following example:

 $ bamCoverage -b reads.bam -o coverage.bw

If you find deepTools useful for your research please cite as:

Ramírez, Fidel, Devon P. Ryan, Björn Grüning, Vivek Bhardwaj, Fabian Kilpert,
Andreas S. Richter, Steffen Heyne, Friederike Dündar,
and Thomas Manke. 2016. “deepTools2: A next Generation Web Server for Deep-Sequencing
Data Analysis.” Nucleic Acids Research, April. doi:10.1093/nar/gkw257.



[ Tools for BAM and bigWig file processing ]
    multiBamSummary         compute read coverages over bam files. Output used for plotCorrelation or plotPCA
    multiBigwigSummary      extract scores from bigwig files. Output used for plotCorrelation or plotPCA
    correctGCBias           corrects GC bias from bam file. Don't use it with ChIP data
    bamCoverage             computes read coverage per bins or regions
    bamCompare              computes log2 ratio and other operations of read coverage of two samples per bins or regions
    bigwigCompare           computes log2 ratio and other operations from bigwig scores of two samples per bins or regions
    computeMatrix           prepares the data from bigwig scores for plotting with plotHeatmap or plotProfile


[ Tools for QC ]
    plotCorrelation         plots heatmaps or scatterplots of data correlation
    plotPCA                 plots PCA
    plotFingerprint         plots the distribution of enriched regions
    bamPEFragmentSize       returns the read length and paired-end distance from a bam file
    computeGCBias           computes and plots the GC bias of a sample
    plotCoverage            plots a histogram of read coverage


[Heatmaps and summary plots]
    plotHeatmap             plots one or multiple heatmaps of user selected regions over different genomic scores
    plotProfile             plots the average profile of user selected regions over different genomic scores
    plotEnrichment          plots the read/fragment coverage of one or more sets of regions


For more information visit: http://deeptools.readthedocs.org
""")

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    return parser


def process_args(args=None):
    args = parse_arguments().parse_args(args)

    return args


def main(args=None):
    if args is None and len(sys.argv) == 1:
        args = ["--help"]
    process_args(args)
