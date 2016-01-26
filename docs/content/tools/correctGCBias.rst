correctGCBias
=============

.. hint:: For background information about the GC bias assessment and correction, see :doc:`computeGCBias`.

.. argparse::
   :ref: deeptools.correctGCBias.parse_arguments
   :prog: correctGCBias

   
Usage example
~~~~~~~~~~~~~~

.. note:: ``correctGCBias`` requires the output of ``computeGCBias``.

.. code:: bash
	
   $ correctGCBias -b H3K27Me3.bam  
      --effectiveGenomeSize 2695000000 
      --genome genome.2bit  
      --GCbiasFrequenciesFile freq_test.txt # output of computeGCBias
      -o gc_corrected.bam

Example output plot
~~~~~~~~~~~~~~~~~~~~

The example shows the GC-bias of a corrected BAM file (output from ``computeGCBias``). 

.. image:: ../../images/test_plots/ExampleCorrectGCBias.png

Galaxy
------

`correctGCBias` is also available in `deepTools Galaxy`_:

.. image:: ../../images/correctGCBias_Galaxy.png 

.. _deepTools Galaxy: http://deeptools2.ie-freiburg.mpg.de/
