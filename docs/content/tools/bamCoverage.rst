bamCoverage
===========

.. argparse::
   :ref: deeptools.bamCoverage.parseArguments
   :prog: bamCoverage


Usage hints
-----------

Command line example using additional options (smaller bin size for higher resolution, normalizing coverage to 1x mouse genome size, excluding chromosome X during the normalization step, and extending reads to an assay specific fragment size of 200 nt):

::

   $ bamCoverage --bam corrected_counts.bam -o Coverage.GCcorrected.SeqDepthNorm.bw
      --binSize 10
      --normalizeTo1x 2150570000
      --ignoreForNormalization chrX
      --fragmentLength 200

* A smaller bin size value will result in a higher resolution of the coverage track but also in a larger file size.
* The 1x normalization (RPGC) requires the input of a value for the **effective genome size**, which is the mappable part of the reference genome. Of course, this value is species specific. The command line help of this tool offers suggestions for a number of model species.
* It might be useful for some studies to exclude certain chromosomes in order to avoid biases, e.g. chromosome X, as male mice contain a pair of each autosome, but usually only a single X chromosome.
* By default, the actual fragment length is estimated from the coordinates of read pairs. The user provided fragment length (e.g. 200 bp) is only used as a fall back for singletons.


Important notes
---------------

* ``--ignoreDuplicates`` : If you already normalized for GC bias using `correctGCbias`, you should absolutely **NOT** set this parameter here!


Galaxy
------

`bamCoverage` is also available in `deepTools Galaxy`_:

.. image:: ../../images/norm_bamCoverage.png 

.. _deepTools Galaxy: http://deeptools.ie-freiburg.mpg.de/