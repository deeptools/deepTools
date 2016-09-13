correctGCBias
=============

.. hint:: For background information about the GC bias assessment and correction, see :doc:`computeGCBias`.

.. argparse::
   :ref: deeptools.correctGCBias.parse_arguments
   :prog: correctGCBias
   :nodefault:

   
Usage example
^^^^^^^^^^^^^^

.. note:: ``correctGCBias`` requires the output of ``computeGCBias`` and a genome file in 2bit format. Most genomes can be found here: http://hgdownload.cse.ucsc.edu/gbdb/. Search for the ``.2bit`` ending. Otherwise, FASTA files can be converted to 2bit using ``faToTwoBit``, which is available here: http://hgdownload.cse.ucsc.edu/admin/exe/

.. code:: bash
	
   $ correctGCBias -b H3K27Me3.bam  
      --effectiveGenomeSize 2695000000 
      --genome genome.2bit  
      --GCbiasFrequenciesFile freq_test.txt # output of computeGCBias
      -o gc_corrected.bam

.. warning:: The GC-corrected BAM file will most likely contain several duplicated reads in regions where the coverage had to increased in order to match the expected read density. This means that you should absolutely avoid using any filtering of duplicate reads during your downstream analyses!
