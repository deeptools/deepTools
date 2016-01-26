plotHeatmap
===========

.. argparse::
   :ref: deeptools.plotHeatmap.parse_arguments
   :prog: plotHeatmap


Usage Example:
~~~~~~~~~~~~~~

The following example creates a heatmap over hg19 transcripts for our test ENCODE datasets. Note that the matrix contains multiple groups of regions (in this case, one for each chromosome used).

.. code:: bash

    plotHeatmap -m matrix_two_groups.gz \
        -out ExampleHeatmap1.png \
        --plotTitle "Test data with default settings"

.. image:: ../../images/test_plots/ExampleHeatmap1.png

plotHeatmap has many options, including the ability to do k-means clustering and change the color map.

.. code:: bash

    plotHeatmap -m matrix_two_groups.gz \
        -out ExampleHeatmap2.png \
        --colorMap Spectral \
        --kmeans 5 \
        --plotTitle "Test data with k-means clustering"

.. image:: ../../images/test_plots/ExampleHeatmap2.png
