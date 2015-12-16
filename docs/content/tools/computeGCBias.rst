computeGCBias
=============

.. argparse::
   :ref: deeptools.computeGCBias.parse_arguments
   :prog: computeGCBias

   
Usage Example:
~~~~~~~~~~~~~~

computeGCBias -b H3K27Me3.bam --effectiveGenomeSize 2695000000 --genome genome.2bit -l 200 -freq freq_test.txt --region X --biasPlot test.gc.plot --plotFileFormat png
