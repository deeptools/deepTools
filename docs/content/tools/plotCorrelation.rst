plotCorrelation
===============

.. contents:: 
    :local:

.. argparse::
   :ref: deeptools.plotCorrelation.parse_arguments
   :prog: plotCorrelation


Example
~~~~~~~~~~~~~~

In the following example, a correlation analysis is performed based on the coverage file computed by :doc:`multiBamCoverage` or :doc:`multiBigwigSummary` for our test ENCODE ChIP-Seq datasets.

**Scatterplot**

Here we make a scatterplot and include pearson correlation coefficients based on the average scores per transcript that were calculated using :doc:`multiBigwigSummary`.

.. code:: bash

    $ deepTools2.0/bin/plotCorrelation \
    -in scores_per_transcript.npz \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation of Average Scores Per Transcript" \
    --whatToPlot scatterplot \
    -o scatterplot_PearsonCorr_bigwigScores.png   \
    --outFileCorMatrix PearsonCorr_bigwigScores.tab 

.. image:: test_plots/scatterplot_PearsonCorr_bigwigScores.png

.. code:: bash

    $ cat PearsonCorr_bigwigScores.tab 
        'H3K27me3'	'H3K4me1'	'H3K4me3'	'HeK9me3'	'input'
        'H3K27me3'	1.0000	-0.1032	-0.1269	-0.0339	-0.0395
        'H3K4me1'	-0.1032	1.0000	0.3985	-0.1863	0.3328
        'H3K4me3'	-0.1269	0.3985	1.0000	-0.0480	0.2822
        'HeK9me3'	-0.0339	-0.1863	-0.0480	1.0000	-0.0353
        'input'	-0.0395	0.3328	0.2822	-0.0353	1.0000


**Heatmap**

Here we plot a heatmap, this time of the Spearman correlation coefficients or read counts that were calculated using :doc:`multiBamCoverage`. 
The dendrogram indicates which samples' read counts are most similar to each other.

.. code:: bash

    $ deepTools2.0/bin/plotCorrelation \
        -in readCounts.npz \
        --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        -o heatmap_SpearmanCorr_readCounts.png   \
        --outFileCorMatrix SpearmanCorr_readCounts.tab 

.. image:: test_plots/heatmap_SpearmanCorr_readCounts.png


plotCorrelation in Galaxy
~~~~~~~~~~~~~~~~~~~~~~~~~

Below is the screenshot showing how to use plotCorrelation with deepTools Galaxy.


.. image:: ../../images/plotCorrelation_galaxy.png
