DeepTools tools
===============

Here, you will find all the options available for the command line
(almost all of them are also available in Galaxy, perhaps named slightly
different).

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

general principles
^^^^^^^^^^^^^^^^^^

-  output format of plots should be indicted by the file ending, e.g.
   MyPlot.pdf will return a pdf, MyPlot.png a png-file
-  all tools that produce plots can also output the underlying data -
   this can be useful in case you don’t like the deepTools visualization
   as you can then use the data matrices produced by deepTools with your
   favorite plotting module, e.g. R or Excel

Parameters to decrease the run time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  numberOfProcessors
-  region - in case you're testing whether a certain plot works and
   gives you the output you're hoping for, you can speed things up by
   focusing on a certain genome region, e.g. chr4 or chr2:100000200000

filtering BAMs while processing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  ignoreDuplicates
-  minMappingQuality
-  samFlagInclude
-  samFlagExclude

+---------------------+
| **If you use        |
| `bamCoverage <https |
| ://github.com/fidel |
| ram/deepTools/wiki/ |
| Normalizations#wiki |
| -bamCoverage>`__    |
| or                  |
| `bamCompare <https: |
| //github.com/fidelr |
| am/deepTools/wiki/N |
| ormalizations#wiki- |
| bamCompare>`__      |
| on samples that     |
| have many           |
| duplicates or many  |
| reads of low        |
| quality, it might   |
| be better to filter |
| the BAM files       |
| *beforehand* as the |
| filtering by        |
| deepTools is done   |
| *after* the scaling |
| factors are         |
| calculated. If you  |
| know that your      |
| files will be       |
| strongly affected   |
| by the filtering,   |
| this might lead to  |
| non-optimal scaling |
| factors!**          |
+---------------------+

To tell a program to use a certain option (e.g. to ignore duplicate
reads), you will have to give the option name preceded by two hyphens
(e.g. --ignoreDuplicates). In the tables on this page, we try to list:

-  the option name as recognized by the program
-  the kind of value that is sometimes expected after the option name
   (see the annotated figure below)
-  a verbose explanation of what the option actually does

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
   tools/bamCorrelate

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