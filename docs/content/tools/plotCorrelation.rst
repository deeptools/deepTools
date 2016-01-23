plotCorrelation
===============

.. contents:: 
    :local:

.. argparse::
   :ref: deeptools.plotCorrelation.parse_arguments
   :prog: plotCorrelation

Background
^^^^^^^^^^^

``plotCorrelation`` computes the overall similarity between two or more files based on read coverage (or other scores) within genomic regions, which must be calculated using either :doc:`multiBamSummary` or :doc:`multiBigwigSummary`.

The result of the correlation computation is a **table of correlation coefficients** that indicates how "strong" the relationship between two samples is and it will consist of numbers between -1 and 1. (-1 indicates perfect anti-correlation, 1 perfect correlation.) 

.. image:: ../../images/QC_bamCorrelate_intro.png

We offer two different functions for the correlation computation: *Pearson* or *Spearman*.

The *Pearson method* measures the **metric differences** between samples and is therefore influenced by outliers.
The *Spearman method* is based on **rankings**.
If you imagine a race with 3 participants where the winner and runner-up are very close together while the third person broke her leg and comes in way, way after the first two, then Pearson would be strongly influenced by the fact that the third person had a great distance to the first ones while Spearman would only care about the fact that person 1 came in first, person 2 came in second and person 3 got the third rank, the distances between them are ignored.

.. tip:: Pearson is an appropriate measure for data that follows a normal distribution, while Spearman does not make this assumption and is generally less driven by outliers, but with the caveat of also being less sensitive.

Here's an example of RNA-seq data from different human cell lines that we had downloaded from https://genome.ucsc.edu/ENCODE/dataMatrix/encodeDataMatrixHuman.html. 

.. image:: ../../images/QC_bamCorrelate_RNAseq.png

As you can see, both correlation calculations more or less agree on which samples are nearly identical (the replicates, indicated by 1 or 2 at the end of the label). The Spearman correlation, however, seems to be more robust and meets our expectations more closely as the two different cell types (HUVEC and IMR90) are clearly separated.

Examples
^^^^^^^^^^^

In the following example, a correlation analysis is performed based on the coverage file computed by :doc:`multiBamSummary` or :doc:`multiBigwigSummary` for our test ENCODE ChIP-Seq datasets.

**Scatterplot**

Here we make pairwose scatterplots of the average scores per transcript that we calculated using :doc:`multiBigwigSummary` and include the Pearson correlation coefficients for each comparison.

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

In addition to scatterplots, heatmaps can be generated where the pairwise correlation coefficients are depicted by varying color intensities and are clustered using hierarchical clustering.

The example here calculates the Spearman correlation coefficients of read counts.
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
^^^^^^^^^^^^^^^^^^^^^^^^^^

Below is the screenshot showing how to use plotCorrelation with deepTools Galaxy.


.. image:: ../../images/plotCorrelation_galaxy.png
