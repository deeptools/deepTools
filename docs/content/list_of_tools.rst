The tools
==================

+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| tool          | type          | input files                       | main output file(s)                    | application                                                                  |
+===============+===============+===================================+========================================+==============================================================================+
| *Correlate    | QC            | 2 or more BAM or bigWig           | table of values                        | Pearson or Spearman correlation between read distributions                   |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
|plotCorrelation| visualization | bam|bigWigCorrelate output        | clustered heatmap                      | visualize the Pearson/Spearman correlation                                   |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| bamFingerprint| QC            | 2 BAM                             | 1 diagnostic plot                      | assess enrichment strength of a ChIP sample                                  |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| computeGCbias | QC            | 1 BAM                             | 2 diagnostic plots                     | calculate the exp. and obs. GC distribution of reads                         |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| correctGCbias | QC            | 1 BAM, output from computeGCbias  | 1 GC-corrected BAM                     | obtain a BAM file with reads distributed according to the genome’s GC conten |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| bamCoverage   | normalization | BAM                               | bedGraph or bigWig                     | obtain the normalized read coverage of a single BAM file                     |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| *Compare      | normalization | 2 BAM or 2 bigWig                 | bedGraph or bigWig                     | normalize 2 files to each other (e.g. log2ratio, difference)                 |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| computeMatrix | visualization | 1 bigWig, 1 BED                   | zipped file for heatmapper or profiler | compute the values needed for heatmaps and summary plots                     |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| heatmapper    | visualization | computeMatrix output              | heatmap of read coverages              | visualize the read coverages for genomic regions                             |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| profiler      | visualization | computeMatrix output              | summary plot (“meta-profile”)          | visualize the average read coverages over a group of genomic regions         |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+

Tools for quality control
^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   tools/bamFingerPrint
   tools/bamCorrelate
   tools/bigwigCorrelate
   tools/computeGCBias

Tools to create tracks
^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   tools/bamCoverage
   tools/bamCompare
   tools/bigwigCompare

Tools for visualization
^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   tools/computeMatrix
   tools/heatmapper
   tools/profiler
   tools/plotCorrelation
   tools/plotCoverage
   tools/plotPCA


General principles
^^^^^^^^^^^^^^^^^^

A typical deepTools command could look like this:

::

    $ /deepTools/bin/bamCoverage --bam myAlignedReads.bam \
    --outFileName myCoverageFile.bigWig \
    --outFileFormat bigwig \
    --fragmentLength 200 \
    --ignoreDuplicates \
    --scaleFactor 0.5

You can always see all available command-line options via --help:

::

    $ /deepTools/bin/bamCoverage --help

-  Output format of plots should be indicated by the file ending, e.g.
   ``MyPlot.pdf`` will return a pdf file, ``MyPlot.png`` a png-file
-  All tools that produce plots can also output the underlying data -
   this can be useful in cases where you don't like the deepTools visualization
   as you can then use the data matrices produced by deepTools with your
   favorite plotting module, e.g. R or Excel
-  The vast majority of command line options are also available in
   Galaxy (in a few cases with minor updates to their naming where needed).

Parameters to decrease the run time
""""""""""""""""""""""""""""""""""""

-  ``numberOfProcessors`` - Number of processors to be used
                        For example, setting ``--numberOfProcessors 10`` will split up the
                        workload internally into 10 chunks, which will be
                        processed in parallel.
-  ``region`` - Allows you to limit the program to a small region.
                        This is particularly useful when you're still trying
                        to figure out the best parameter setting, e.g., for 
                        certain plots. You can focus on a certain genome
                        region by setting, e.g., ``--region chr2`` or even
                        ``--region chr2:100000-200000``

These parameters are optional and available throughout almost all deepTools.

Filtering BAMs while processing
""""""""""""""""""""""""""""""""

-  ``ignoreDuplicates`` 
                        Reads with the same orientation and start
                        position will be considered only once. If reads are
                        paired, the mate is also evaluated
-  ``minMappingQuality``
                        Only reads with a mapping quality score equal
                        or higher than the specified value are considered
-  ``samFlagInclude``
                        Include reads based on the SAM flag, e.g.
                        ``--samFlagInclude 64`` gets reads that are first in
                        a pair. For translating SAM flags into English, go to:
                        https://broadinstitute.github.io/picard/explain-flags.html
-  ``samFlagExclude``
                        Exclude reads based on the SAM flags - see previous explanation.

These parameters are optional and available throughout deepTools.

.. warning::  If you know that your files will be strongly affected by the filtering
 of duplicates or reads of low quality, you should consider removing
 those reads *before* using bamCoverage or bamCompare as the filtering
 by deepTools is done *after* the scaling factors are calculated!