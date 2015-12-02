deepTools overview
==================

Here, you will find all command line options. They are also available in
Galaxy (in a few cases with minor updates to their naming where needed).

You can always see all available command-line options via --help:

::

    $ /deepTools/bin/bamCoverage --help

A typical deepTools command could look like this:

::

    $ /deepTools/bin/bamCoverage --bam myAlignedReads.bam \
    --outFileName myCoverageFile.bigWig \
    --outFileFormat bigwig \
    --fragmentLength 200 \
    --ignoreDuplicates \
    --scaleFactor 0.5

General principles
^^^^^^^^^^^^^^^^^^

-  Output format of plots should be indicated by the file ending, e.g.
   MyPlot.pdf will return a pdf, MyPlot.png a png-file
-  All tools that produce plots can also output the underlying data -
   this can be useful in cases where you don't like the deepTools visualization
   as you can then use the data matrices produced by deepTools with your
   favorite plotting module, e.g. R or Excel

Parameters to decrease the run time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  ``numberOfProcessors`` - Number of processors to be used, e.g.
                        ``--numberOfProcessors 10`` will split up the
                        workload internally into 10 chunks which will be
                        processed in parallel
-  ``region`` - In case you're testing whether a certain plot works and
                        gives you the output you're hoping for, you can
                        speed things up by focusing on a certain genome
                        region, e.g. ``--region chr2`` or even
                        ``--region chr2:100000-200000``

These parameters are optional and available throughout deepTools.

Filtering BAMs while processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  ``ignoreDuplicates`` - Reads with the same orientation and start
                        position will be considered only once. If reads are
                        paired, the mate is also evaluated
-  ``minMappingQuality`` - Only reads with a mapping quality score equal
                        or higher than the specified value are considered
-  ``samFlagInclude`` - Include reads based on the SAM flag, e.g.
                        ``--samFlagInclude 64`` gets reads that are first in
                        pair. A good place helping to unravel SAM flags is
                        https://broadinstitute.github.io/picard/explain-flags.html
-  ``samFlagExclude`` - Exclude reads based on the SAM flags

These parameters are optional and available throughout deepTools.


.. warning::  If you know that your files will be strongly affected by the filtering
 of duplicates or reads of low quality, you should consider removing
 those reads *before* using bamCoverage or bamCompare as the filtering
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

Tools for quality control
^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   tools/bamFingerPrint
   tools/bamCorrelate
   tools/computeGCBias

Tools to create tracks
^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   tools/bamCoverage
   tools/bamCompare

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

