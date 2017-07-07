estimateReadFiltering
=====================

.. contents:: 
    :local:

.. argparse::
   :ref: deeptools.estimateReadFiltering.parseArguments
   :prog: estimateReadFiltering
   :nodefault:


Background
^^^^^^^^^^

Many tools within deepTools allow one to filter BAM files according to alignment mapping qualities or other criteria. It's difficult to know ahead of time how these various settings will affect the number of filtered reads. Consequently, ``estimateReadFiltering`` can be used to approximate the number of reads in a BAM file or files that will be filtered according to one or more criteria. This can also be used the quickly estimate the duplication level in a BAM file.

Usage example
^^^^^^^^^^^^^

``estimateReadFiltering`` needs one or more sorted and indexed BAM files and the desired filtering criteria.

.. code:: bash

    $ estimateReadFiltering -b paired_chr2L.bam \
    --minMappingQuality 5 --samFlagInclude 16 \
    --samFlagExclude 256 --ignoreDuplicates

By default, the output is printed to the screen. You can change this with the ``-o`` option. The output is a tab-separated file:

    Sample  Total Reads     Mapped Reads    Alignments in blacklisted regions       Estimated mapped reads filtered Below MAPQ      Missing Flags   Excluded Flags  Internally-determined Duplicates        Marked Duplicates  Singletons      Wrong strand
    paired_chr2L.bam        12644   12589   0       6313.2  4114.0  6340.0  0.0     1163.0  0.0     55.0    0.0

The columns are as follows:

 * Total reads (including unmapped)
 * Unmapped reads
 * Reads in blacklisted regions (--blackListFileName)

The following metrics are estimated according to the --binSize and --distanceBetweenBins parameters
 * Estimated mapped reads filtered (the total number of mapped reads filtered for any reason)
 * Alignments with a below threshold MAPQ (--minMappingQuality)
 * Alignments with at least one missing flag (--samFlagInclude)
 * Alignments with undesirable flags (--samFlagExclude)
 * Duplicates determined by deepTools (--ignoreDuplicates)
 * Duplicates marked externally (e.g., by picard)
 * Singletons (paired-end reads with only one mate aligning)
 * Wrong strand (due to --filterRNAstrand)

The sum of these may be more than the total number of reads. Note that alignments are sampled from bins of size --binSize spaced --distanceBetweenBins apart.

