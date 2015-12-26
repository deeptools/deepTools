bamFingerprint
==========================

.. argparse::
   :ref: deeptools.bamFingerprint.parse_arguments
   :prog: bamFingerprint


Usage Example:
~~~~~~~~~~~~~~

The following example plots the the fingerprint for our test ENCODE
Chip-Seq datasets, for chromosome X. Here we are skipping the regions
with zero reads using the *--skipZeros* option.

.. code:: bash

    bamFingerprint -p 25 -b testDataset/H3K27Me3.bam testDataset/H3K4Me1.bam testDataset/H3K4Me3.bam testDataset/H3K9Me3.bam testDataset/Input.bam testDataset/RNAseq.bam -r X --skipZeros -plot test_plots/bamFP_chrX.png

.. image:: test_plots/bamFP_chrX.png


bamFingerprint on Galaxy:
~~~~~~~~~~~~~~~~~~~~~~~~~

Below is the screenshot showing how to use bamFingerprint on the deeptools galaxy.


.. image:: ../../images/bamFP_galaxy.png
