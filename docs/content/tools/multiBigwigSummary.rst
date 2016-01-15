bigwigCorrelate
===============

.. argparse::
   :ref: deeptools.multiBigwigSummary.parse_arguments
   :prog: bigwigCorrelate


Usage Example:
~~~~~~~~~~~~~~

The following example computes the correlation for our test ENCODE
Chip-Seq datasets.

.. code:: bash

   multiBigwigSummary bins -p 20 -b testDataset/H3K4Me1.bigWig \
      testDataset/H3K4Me3.bigWig testDataset/H3K27Me3.bigWig \
      testDataset/Input.bigWig -out histoneMarks_bigwig_corr.npz

The output matrix file can now be plotted using :doc:`plotCorrelation`.


bigwigCorrelate on Galaxy:
~~~~~~~~~~~~~~~~~~~~~~~~~

Below is the screenshot showing how to use multiBigwigSummary on the deeptools galaxy.


.. image:: ../../images/bigwiCorr_galaxy.png
