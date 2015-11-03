bamCoverage
===========

.. argparse::
   :ref: deeptools.bamCoverage.parseArguments
   :prog: bamCoverage

Usage hints
-----------

::

	$ bamCoverage --bam corrected_counts.bam --binSize 10 \
		--normalizeTo1x 2150570000 --fragmentLength 200 \
		-o Coverage.GCcorrected.SeqDepthNorm.bw --ignoreForNormalization chrX

* The bin size can be chosen completely to your liking (`-bs` option). The smaller it is, the bigger the resulting file will be.
* The above shown example was for a mouse sample, therefore the effective genome size for mouse had to be indicated once it was decided that the file should be normalize to 1x coverage.
* Chromosome X was excluded from sampling the regions for normalization as the sample was from a male mouse that therefore contained pairs of autosomes, but only a single X chromosome.
* The fragment length of 200 bp is only the fall-back option of `bamCoverage` as the sample provided here was done with paired-end sequencing. `bamCoverage` will resort to the user-specified fragment length only if it encounters singletons.
* `--ignoreDuplicates` - important note! if you normalized for GC bias using `correctGCbias`, you should absolutely **NOT** set this parameter 
 
Here's an example of how `bamCoverage` can be used with `deepTools Galaxy`_:

.. image:: ../../images/norm_bamCoverage.png 

.. _deepTools Galaxy: http://deeptools.ie-freiburg.mpg.de/