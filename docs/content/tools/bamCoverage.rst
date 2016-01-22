bamCoverage
===========

.. contents:: 
    :local:

.. argparse::
   :ref: deeptools.bamCoverage.parseArguments
   :prog: bamCoverage


Usage hints
-----------

* A smaller bin size value will result in a higher resolution of the coverage track but also in a larger file size.
* The 1x normalization (RPGC) requires the input of a value for the **effective genome size**, which is the mappable part of the reference genome. Of course, this value is species-specific. The command line help of this tool offers suggestions for a number of model species.
* It might be useful for some studies to exclude certain chromosomes in order to avoid biases, e.g. chromosome X, as male mice contain a pair of each autosome, but usually only a single X chromosome.
* By default, the read length is **NOT** extended! This is the preferred setting for **spliced-read** data like RNA-seq, where one usually wants to rely on the detected read locations only. A read extension would neglect potential splice sites in the unmapped part of the fragment.
  Other data, e.g. Chip-seq, where fragments are known to map contiguously, should be processed with read extension (``--extendReads [INT]``).
* For paired-end data, the fragment length is generally defined by the two read mates. The user provided fragment length is only used as a fallback for singletons or mate reads that map too far apart (with a distance greater than four times the fragment length or are located on different chromosomes).

.. warning:: If you already normalized for GC bias using `correctGCbias`, you should absolutely **NOT** set the parameter ``--ignoreDuplicates``!


Usage examples for ChIP-seq
---------------------------

This is an example for ChIP-seq data using additional options (smaller bin size for higher resolution, normalizing coverage to 1x mouse genome size, excluding chromosome X during the normalization step, and extending reads):

.. code:: bash

    bamCoverage --bam a.bam -o a.SeqDepthNorm.bw \
        --binSize 10
        --normalizeTo1x 2150570000
        --ignoreForNormalization chrX
        --extendReads


`bamCoverage` is also available in `deepTools Galaxy`_:

.. image:: ../../images/norm_bamCoverage.png 

.. _deepTools Galaxy: http://deeptools.ie-freiburg.mpg.de/


Usage examples for RNA-seq
--------------------------

Note that some BAM files are filtered based on SAM flags (`Explain SAM flags <https://broadinstitute.github.io/picard/explain-flags.html>`_).


.. code:: bash

    # Regular bigWig track
    bamCoverage -b a.bam -o a.bw

    # Forward strand only (for paired-end stranded library)
    samtools view -b -f 128 -F 16 a.bam > a.fwd1.bam
    samtools view -b -f 64 -F 32 a.bam > a.fwd2.bam
    samtools merge -f fwd.bam fwd1.bam fwd2.bam
    bamCoverage -b fwd.bam -o a.fwd.bw
    rm *.fwd*.bam

    # Reverse strand only (for paired-end stranded library)
    samtools view -b -f 144 a.bam > a.rev1.bam
    samtools view -b -f 96 a.bam > a.rev2.bam
    samtools merge -f rev.bam rev1.bam rev2.bam
    bamCoverage -b rev.bam -o a.rev.bw
    rm *.rev*.bam

    # Forward strand only (for single-end stranded library)
    bamCoverage -b a.bam -o a.fwd.bw --samFlagExclude 16

    # Reverse strand only (for single-end stranded library)
    bamCoverage -b a.bam -o a.rev.bw --samFlagInclude 16