bamFilterReads
==============

.. contents:: 
    :local:

.. argparse::
   :ref: deeptools.bamFilterReads.parseArguments
   :prog: bamFilterReads
   :nodefault:


Background
^^^^^^^^^^

This tool estimates the number of alignments that would be excluded from one or more BAM files given a variety of filtering criteria. This is useful for estimating the duplication rate in an experiment or more generally seeing what the effect of various option choices will be in other deepTools tools without actually spending the time to run them.

Usage example
^^^^^^^^^^^^^

``bamFilterReads`` needs a sorted and indexed BAM file and the desired filtering criteria.

.. code:: bash

    $ bamFilterReads -b paired_chr2L.bam \
    --minMappingQuality 5 --samFlagInclude 16 \
    --samFlagExclude 256 --ignoreDuplicates \
    -o filtered.bam --filterMetrics metrics.txt

The filtered results are then written to the file specified by ``-o``. If you would like to store metrics about the number of reads seen and the number remaining after filtering, then specify the file for that with ``--filterMetrics``. An example metrics file is below:

    #bamFilterReads --filterMetrics
    #File	Reads Remaining	Total Initial Reads
    paired_chr2L.bam	8440	12644
