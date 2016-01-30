plotHeatmap
===========

.. contents:: 
    :local:

.. argparse::
   :ref: deeptools.plotHeatmap.parse_arguments
   :prog: plotHeatmap


Details
^^^^^^^^



The following tables summarizes the kinds of optional outputs that are available with the three tools.

+-----------------------------------+--------------------------------+-------------------+-----------------+-----------------+
|  **optional output type**         | **command**                    | **computeMatrix** | **plotHeatmap** | **plotProfile** |
+-----------------------------------+--------------------------------+-------------------+-----------------+-----------------+
| values underlying the heatmap     | ``--outFileNameMatrix``        | yes               | yes             | no              |
+-----------------------------------+--------------------------------+-------------------+-----------------+-----------------+
| values underlying the profile     | ``--outFileNameData``          | no                | yes             | yes             |
+-----------------------------------+--------------------------------+-------------------+-----------------+-----------------+
| sorted and/or filtered regions    | ``--outFileSortedRegions``     | yes               | yes             | yes             |
+-----------------------------------+--------------------------------+-------------------+-----------------+-----------------+


Usage examples
^^^^^^^^^^^^^^

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
