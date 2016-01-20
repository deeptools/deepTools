The tools
=========

.. note:: With the release of deepTools 2.0, we renamed a couple of tools: **heatmapper** to :doc:`tools/plotHeatmap`; **profiler** to :doc:`tools/plotProfile`; **bamCorrelate** to :doc:`tools/multiBamCoverage`; **bigwigCorrelate** to :doc:`tools/multiBigwigSummary` and **bamFingerprint** to :doc:`tools/plotFingerprint`. For more, see :doc:`changelog`.

+--------------------------------+------------------+-----------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
| tool                           | type             | input files                       | main output file(s)                        | application                                                                       |
+================================+==================+===================================+============================================+===================================================================================+
|:doc:`tools/multiBamCoverage`   | data integration | 2 or more BAM                     | interval-based table of values             | perform cross-sample analyses of read counts --> plotCorrelation, plotPCA         |
+--------------------------------+------------------+-----------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/multiBigwigSummary` | data integration | 2 or more bigWig                  | interval-based table of values             |  perform cross-sample analyses of genome-wide scores --> plotCorrelation, plotPCA |
+--------------------------------+------------------+-----------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/plotCorrelation`    | visualization    | bam/multiBigwigSummary output     | clustered heatmap                          | visualize the Pearson/Spearman correlation                                        |
+--------------------------------+------------------+-----------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/plotPCA`            | visualization    | bam/multiBigwigSummary output     | 2 PCA plots                                | visualize the principal component analysis                                        |
+--------------------------------+------------------+-----------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/plotFingerprint`    | QC               | 2 BAM                             | 1 diagnostic plot                          | assess enrichment strength of a ChIP sample                                       |
+--------------------------------+------------------+-----------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/computeGCBias`      | QC               | 1 BAM                             | 2 diagnostic plots                         | calculate the exp. and obs. GC distribution of reads                              |
+--------------------------------+------------------+-----------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/correctGCBias`      | QC               | 1 BAM, output from computeGCbias  | 1 GC-corrected BAM                         | obtain a BAM file with reads distributed according to the genome’s GC content     |
+--------------------------------+------------------+-----------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/bamCoverage`        | normalization    | BAM                               | bedGraph or bigWig                         | obtain the normalized read coverage of a single BAM file                          |
+--------------------------------+------------------+-----------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/bamCompare`         | normalization    | 2 BAM                             | bedGraph or bigWig                         | normalize 2 files to each other (e.g. log2ratio, difference)                      |
+--------------------------------+------------------+-----------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/computeMatrix`      | data integration | 1 or more bigWig, 1 or more BED   | zipped file for plotHeatmap or plotProfile | compute the values needed for heatmaps and summary plots                          |
+--------------------------------+------------------+-----------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/plotHeatmap`        | visualization    | computeMatrix output              | heatmap of read coverages                  | visualize the read coverages for genomic regions                                  |
+--------------------------------+------------------+-----------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/plotProfile`        | visualization    | computeMatrix output              | summary plot (“meta-profile”)              | visualize the average read coverages over a group of genomic regions              |
+--------------------------------+------------------+-----------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/plotCoverage`       | visualization    | 1 or more bam                     | 2 diagnostic plots                         | visualize the average read coverages over sampled genomic  positions              |
+--------------------------------+------------------+-----------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/bamPEFragmentSize`  | information      | 1  BAM                            | text with paired-end fragment length       | obtain the average fragment length from paired ends                               |
+--------------------------------+------------------+-----------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+

.. contents:: 
    :local:

General principles
^^^^^^^^^^^^^^^^^^

A typical deepTools command could look like this:

.. code:: bash

    $ bamCoverage --bam myAlignedReads.bam \
    --outFileName myCoverageFile.bigWig \
    --outFileFormat bigwig \
    --fragmentLength 200 \
    --ignoreDuplicates \
    --scaleFactor 0.5

You can always see all available command-line options via --help:

.. code:: bash

    $ bamCoverage --help

-  Output format of plots should be indicated by the file ending, e.g.
   ``MyPlot.pdf`` will return a pdf file, ``MyPlot.png`` a png-file
-  All tools that produce plots can also output the underlying data -
   this can be useful in cases where you don't like the deepTools visualization,
   as you can then use the data matrices produced by deepTools with your
   favorite plotting tool, such as R
-  The vast majority of command line options are also available in
   Galaxy (in a few cases with minor changes to their naming).

Parameters to decrease the run time
"""""""""""""""""""""""""""""""""""

-  ``numberOfProcessors`` - Number of processors to be used
                        For example, setting ``--numberOfProcessors 10`` will split up the
                        workload internally into 10 chunks, which will be
                        processed in parallel.
-  ``region`` - Process only a single genomic region.
                        This is particularly useful when you're still trying
                        to figure out the best parameter setting, e.g., for 
                        certain plots. You can focus on a certain genomic
                        region by setting, e.g., ``--region chr2`` or 
                        ``--region chr2:100000-200000``

These parameters are optional and available throughout almost all deepTools.

Filtering BAMs while processing
"""""""""""""""""""""""""""""""

-  ``ignoreDuplicates`` 
                        Reads with the same orientation and start
                        position will be considered only once. If reads are
                        paired, the mate is also evaluated
-  ``minMappingQuality``
                        Only reads with a mapping quality score of at least this are considered
-  ``samFlagInclude``
                        Include reads based on the SAM flag, e.g.
                        ``--samFlagInclude 64`` gets reads that are first in
                        a pair. For translating SAM flags into English, go to:
                        https://broadinstitute.github.io/picard/explain-flags.html
-  ``samFlagExclude``
                        Exclude reads based on the SAM flags - see previous explanation.

These parameters are optional and available throughout deepTools.

.. warning::  If you know that your files will be strongly affected by the filtering
 of duplicates or reads of low quality then consider removing
 those reads *before* using bamCoverage or bamCompare, as the filtering
 by deepTools is done *after* the scaling factors are calculated!

On the command line, to tell a program to use a certain option
(e.g. to ignore duplicate reads), you will have to give the option name
preceded by two hyphens (e.g. ``--ignoreDuplicates``).

The tables on this page list:

-  The option name as recognized by the program
-  The kind of value that is sometimes expected after the option name
   (see the annotated figure below)
-  A verbose explanation of what the option actually does

The texts here are adjusted for readability, they might not match the
help text that you see in the command line word by word.


Tools for BAM and bigWig file processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:doc:`tools/multiBamCoverage`
"""""""""""""""""""""""""""""""
:doc:`tools/multiBigwigSummary`
"""""""""""""""""""""""""""""""
:doc:`tools/correctGCBias`
"""""""""""""""""""""""""""""""
:doc:`tools/bamCoverage`
"""""""""""""""""""""""""""""""
:doc:`tools/bamCompare`
"""""""""""""""""""""""""""""""
:doc:`tools/bigwigCompare`
"""""""""""""""""""""""""""""""
:doc:`tools/computeMatrix`
"""""""""""""""""""""""""""""""

Tools for QC
^^^^^^^^^^^^^

:doc:`tools/plotCorrelation`
"""""""""""""""""""""""""""""""
:doc:`tools/plotPCA`
"""""""""""""""""""""""""""""""
:doc:`tools/plotFingerprint`
"""""""""""""""""""""""""""""""
:doc:`tools/bamPEFragmentSize`
"""""""""""""""""""""""""""""""
:doc:`tools/computeGCBias`
"""""""""""""""""""""""""""""""
:doc:`tools/plotCoverage`
"""""""""""""""""""""""""""""""

Heatmaps and summary plots
^^^^^^^^^^^^^^^^^^^^^^^^^^

:doc:`tools/plotHeatmap`
"""""""""""""""""""""""""""""""
:doc:`tools/plotProfile`
"""""""""""""""""""""""""""""""
