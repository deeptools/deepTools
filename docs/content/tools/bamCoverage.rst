bamCoverage
===========

.. argparse::
   :ref: deeptools.bamCoverage.parseArguments
   :prog: bamCoverage


Usage hints
-----------

This is an example using additional options (smaller bin size for higher resolution, normalizing coverage to 1x mouse genome size, excluding chromosome X during the normalization step, and extending reads):

::

   $ bamCoverage --bam corrected_counts.bam -o Coverage.GCcorrected.SeqDepthNorm.bw
      --binSize 10
      --normalizeTo1x 2150570000
      --ignoreForNormalization chrX
      --extendReads

* A smaller bin size value will result in a higher resolution of the coverage track but also in a larger file size.
* The 1x normalization (RPGC) requires the input of a value for the **effective genome size**, which is the mappable part of the reference genome. Of course, this value is species-specific. The command line help of this tool offers suggestions for a number of model species.
* It might be useful for some studies to exclude certain chromosomes in order to avoid biases, e.g. chromosome X, as male mice contain a pair of each autosome, but usually only a single X chromosome.
* By default, the read length is **NOT** extended! This is the preferred setting for **spliced-read** data like RNA-seq, where one usually wants to rely on the detected read locations only. A read extension would neglect potential splice sites in the unmapped part of the fragment.
  Other data, e.g. Chip-seq, where fragments are known to map contiguously, should be processed with read extension (``--extendReads [INT]``).
* For paired-end data, the fragment length is generally defined by the two read mates. The user provided fragment length is only used as a fallback for singletons or mate reads that map too far apart (with a distance greater than four times the fragment length or are located on different chromosomes).


Important notes
---------------

* ``--ignoreDuplicates`` : If you already normalized for GC bias using `correctGCbias`, you should absolutely **NOT** set this parameter here!


Galaxy
------

`bamCoverage` is also available in `deepTools Galaxy`_:

.. image:: ../../images/norm_bamCoverage.png 

.. _deepTools Galaxy: http://deeptools.ie-freiburg.mpg.de/
