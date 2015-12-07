bigwigCorrelate
==========================

.. argparse::
   :ref: deeptools.bigwigCorrelate.parse_arguments
   :prog: bigwigCorrelate


Usage Example:
~~~~~~~~~~~~~~

The following example computes the correlation for our test ENCODE
Chip-Seq datasets.

.. code:: bash

    bigwigCorrelate bins -p 20 -b testDataset/H3K4Me1.bigWig testDataset/H3K4Me3.bigWig testDataset/H3K27Me3.bigWig testDataset/Input.bigWig -out histoneMarks_bigwig_corr.npz

The output matrix file can now be plotted using `plotCorrelation <>`__
tool.
