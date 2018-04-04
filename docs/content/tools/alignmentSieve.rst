alignmentSieve
==============

.. contents:: 
    :local:

.. argparse::
   :ref: deeptools.alignmentSieve.parseArguments
   :prog: alignmentSieve
   :nodefault:


Background
^^^^^^^^^^

This tool filters alignments in a BAM/CRAM file according the the specified parameters. It can optionally output to BEDPE format, possibly with the fragment ends shifted in a custom manner.

Usage example
^^^^^^^^^^^^^

``alignmentSieve`` needs a sorted and indexed BAM file and the desired filtering criteria.

.. code:: bash

    $ alignmentSieve -b paired_chr2L.bam \
    --minMappingQuality 5 --samFlagInclude 16 \
    --samFlagExclude 256 --ignoreDuplicates \
    -o filtered.bam --filterMetrics metrics.txt

The alignments passing the filtering criteria are then written to the file specified by ``-o``. You can additionally save alignments **NOT** passing the filtering criteria with the ``-filteredOutReads`` If you would like to store metrics about the number of reads seen and the number remaining after filtering, then specify the file for that with ``--filterMetrics``. An example metrics file is below:

    #bamFilterReads --filterMetrics
    #File	Reads Remaining	Total Initial Reads
    paired_chr2L.bam	8440	12644

Instead of a BAM file, a BEDPE file (suitable for input into MACS2) can be produced. Like the BAM/CRAM output, BEDPE also allows shifting of fragment ends, as is often desirable in ATAC-seq and related protocols:

.. code:: bash

    $ alignmentSieve -b paired_chr2L.bam \
    --minFragmentLength 140 --BED \
    --shift -5 3 -o fragments.bedpe

The ``--shift`` option can take either 2 or 4 integers. If two integers are given, then the first value shifts the left-most end of a fragment and the second the right-most end. Positive values shift to the right and negative values to the left. See below for how the above settings would shift a single fragment::

         ----> read 1
                     read 2 <----

         ------------------------ fragment
    
    -------------------------------- shifted fragment

The same results will be produced if read 1 and read 2 are swapped. If, instead, the protocol is strand-specific, then the first set of integers in a pair would be applied to fragments where read 1 precedes read 2, and the second set to cases where read 2 precedes read 1. In this case, the first value in each pair is applied to the end of read 1 and the second to the end of read 2. Take the following command as an example:

.. code:: bash

    $ alignmentSieve -b paired_chr2L.bam \
    --minFragmentLength 140 --BED \
    --shift -5 3 -1 4 -o fragments.bedpe

Given that, the ``-5 3`` set would produce the following::

         ----> read 1
                     read 2 <----

         ------------------------ fragment
    
    -------------------------------- shifted fragment

and the ``-1 4`` set would produce the following::

         ----> read 2
                     read 1 <----

         ------------------------ fragment

             --------------------- shifted fragment

As can be seen, such fragments are considered to be on the ``-`` strand, so negative values then shift to the left on its frame of reference (thus, to the right relative to the ``+`` strand).

.. note::
    If the ``--shift`` or ``--ATACshift`` options are used, then only properly-paired reads will be used.
