bamCoverage
===========

.. contents::
    :local:

.. image:: ../../images/norm_IGVsnapshot_indFiles.png

If you are not familiar with BAM, bedGraph and bigWig formats, you can read up on that in our :doc:`../help_glossary`

.. argparse::
   :ref: deeptools.bamCoverage.parseArguments
   :prog: bamCoverage
   :nodefault:


Usage hints
^^^^^^^^^^^^

* A smaller bin size value will result in a higher resolution of the coverage track but also in a larger file size.
* The 1x normalization (RPGC) requires the input of a value for the **effective genome size**, which is the mappable part of the reference genome. Of course, this value is species-specific. The command line help of this tool offers suggestions for a number of model species.
* It might be useful for some studies to exclude certain chromosomes in order to avoid biases, e.g. chromosome X, as male mice contain a pair of each autosome, but usually only a single X chromosome.
* By default, the read length is **NOT** extended! This is the preferred setting for **spliced-read** data like RNA-seq, where one usually wants to rely on the detected read locations only. A read extension would neglect potential splice sites in the unmapped part of the fragment.
  Other data, e.g. Chip-seq, where fragments are known to map contiguously, should be processed with read extension (``--extendReads [INTEGER]``).
* For paired-end data, the fragment length is generally defined by the two read mates. The user provided fragment length is only used as a fallback for singletons or mate reads that map too far apart (with a distance greater than four times the fragment length or are located on different chromosomes).

.. warning:: If you already normalized for GC bias using ``correctGCbias``, you should absolutely **NOT** set the parameter ``--ignoreDuplicates``!

.. note:: Like BAM files, bigWig files are compressed, binary files. If you would like to see the coverage values, choose the bedGraph output via ``--outFileFormat``.

Usage example for ChIP-seq
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is an example for ChIP-seq data using additional options (smaller bin size for higher resolution, normalizing coverage to 1x mouse genome size, excluding chromosome X during the normalization step, and extending reads):

.. code:: bash

    bamCoverage --bam a.bam -o a.SeqDepthNorm.bw \
        --binSize 10
        --normalizeUsing RPGC
        --effectiveGenomeSize 2150570000
        --ignoreForNormalization chrX
        --extendReads

If you had run the command with ``--outFileFormat bedgraph``, you could easily peak into the resulting file.

.. code:: bash

    $ head SeqDepthNorm_chr19.bedgraph
    19	60150	60250	9.32
    19	60250	60450	18.65
    19	60450	60650	27.97
    19	60650	60950	37.29
    19	60950	61000	27.97
    19	61000	61050	18.65
    19	61050	61150	27.97
    19	61150	61200	18.65
    19	61200	61300	9.32
    19	61300	61350	18.65

As you can see, each row corresponds to one region. If consecutive bins have the same number of reads overlapping, they will be merged.

Usage examples for RNA-seq
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Note that some BAM files are filtered based on SAM flags (`Explain SAM flags <https://broadinstitute.github.io/picard/explain-flags.html>`_).

Regular bigWig track
~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    bamCoverage -b a.bam -o a.bw


Separate tracks for each strand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes it makes sense to generate two independent :ref:`bigWig` files for all reads on the forward and reverse strand, respectively.
As of deepTools version 2.2, one can simply use the ``--filterRNAstrand`` option, such as ``--filterRNAstrand forward`` or ``--filterRNAstrand reverse``.
This handles paired-end and single-end datasets. For older versions of deepTools, please see the instructions below.

.. note:: The ``--filterRNAstrand`` option assumes the sequencing library generated from ILLUMINA dUTP/NSR/NNSR methods, which are the most commonly used method for
          library preparation, where Read 2 (R2) is in the direction of RNA strand (**reverse-stranded** library). However other methods exist, which generate read
          R1 in the direction of RNA strand (`see this review <http://www.nature.com/nmeth/journal/v7/n9/full/nmeth.1491.html>`_). For these libraries,
          ``--filterRNAstrand`` will have an opposite behavior, i.e. ``--filterRNAstrand forward`` will give you reverse strand signal and vice-versa.

Versions before 2.2
*******************

To follow the examples, you need to know that ``-f`` will tell ``samtools view`` to **include** reads with the indicated flag, while ``-F`` will lead to the **exclusion** of reads with the respective flag.

**For a stranded `single-end` library**

.. code:: bash

    # Forward strand
    bamCoverage -b a.bam -o a.fwd.bw --samFlagExclude 16

    # Reverse strand
    bamCoverage -b a.bam -o a.rev.bw --samFlagInclude 16



**For a stranded `paired-end` library**

Now, this gets a bit cumbersome, but future releases of deepTools will make this more straight-forward.
For now, bear with us and perhaps read up on SAM flags, e.g. `here <http://ppotato.wordpress.com/2010/08/25/samtool-bitwise-flag-paired-reads/>`_.

For paired-end samples, we assume that a proper pair should have the mates on opposing strands where the Illumina strand-specific protocol produces reads in a ``R2-R1`` orientation. We basically follow the recipe given `in this biostars tutorial <https://www.biostars.org/p/92935/>`_.

To get the file for transcripts that originated from the **forward strand**:

.. code:: bash


    # include reads that are 2nd in a pair (128);
    # exclude reads that are mapped to the reverse strand (16)
    $ samtools view -b -f 128 -F 16 a.bam > a.fwd1.bam

    # exclude reads that are mapped to the reverse strand (16) and
    # first in a pair (64): 64 + 16 = 80
    $ samtools view -b -f 80 a.bam > a.fwd2.bam

    # combine the temporary files
    $ samtools merge -f fwd.bam a.fwd1.bam a.fwd2.bam

    # index the filtered BAM file
    $ samtools index fwd.bam

    # run bamCoverage
    $ bamCoverage -b fwd.bam -o a.fwd.bigWig

    # remove the temporary files
    $ rm a.fwd*.bam

To get the file for transcripts that originated from the **reverse strand**:

.. code:: bash

    # include reads that map to the reverse strand (128)
    # and are second in a pair (16): 128 + 16 = 144
    $ samtools view -b -f 144 a.bam > a.rev1.bam

    # include reads that are first in a pair (64), but
    # exclude those ones that map to the reverse strand (16)
    $ samtools view -b -f 64 -F 16 a.bam > a.rev2.bam

    # merge the temporary files
    $ samtools merge -f rev.bam a.rev1.bam a.rev2.bam

    # index the merged, filtered BAM file
    $ samtools index rev.bam

    # run bamCoverage
    $ bamCoverage -b rev.bam -o a.rev.bw

    # remove temporary files
    $ rm a.rev*.bam
