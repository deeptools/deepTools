plotCorrelation
==========================

.. argparse::
   :ref: deeptools.bamFingerprint.parse_arguments
   :prog: bamFingerprint


Usage Example:
~~~~~~~~~~~~~~

The following example plots the correlation matrix computed by
`bamCorrelate <>`__ for our test ENCODE Chip-Seq datasets.

**Scatterplot**

Here we are plotting a scatterplot with pearson correlation.

.. code:: bash

    plotCorrelation -in testDatset-results/histoneMarks_bigwig_corr.npz -o histoneMarks_corr-scatter.png -T "test data correlations" -p scatterplot --removeOutliers -c pearson

**Heatmap**

Here we are plotting a heatmap, this time using spearman correlation.

.. code:: bash

    plotCorrelation -in testDatset-results/histoneMarks_bigwig_corr.npz -o histoneMarks_corr-heatmap.png -T "test data correlations" -p heatmap --removeOutliers -c spearman
