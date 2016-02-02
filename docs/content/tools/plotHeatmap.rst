plotHeatmap
===========

.. contents:: 
    :local:

.. argparse::
   :ref: deeptools.plotHeatmap.parse_arguments
   :prog: plotHeatmap


Details
^^^^^^^^

``plotHeatmap`` does not change the values that ``computeMatrix`` calculated, it simply translates them into heatmaps and summary plots.
It offers a large variety of parameters to explore various visualizations and customize the resulting image (see the commands above).

In addition, you can retrieve all the data tables underlying the various plots including the regions that were used to generate the final plot. 
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

.. tip:: For more details on the optional output, see the examples for :doc:`computeMatrix`.


Usage examples
^^^^^^^^^^^^^^

The following example creates a heatmap over hg19 transcripts for our test ENCODE datasets. Note that the matrix contains multiple groups of regions (in this case, one for each chromosome used).

.. code:: bash

   $ plotHeatmap -m matrix_two_groups.gz \
        -out ExampleHeatmap1.png \
        --plotTitle "Test data with default settings"

.. image:: ../../images/test_plots/ExampleHeatmap1.png

As mentioned above, ``plotHeatmap`` has many options, including the ability to do k-means clustering and change the color map.

.. code:: bash

   $ plotHeatmap -m matrix_two_groups.gz \
        -out ExampleHeatmap2.png \
        --colorMap Spectral \
        --kmeans 5 \
        --plotTitle "Test data with k-means clustering"

.. image:: ../../images/test_plots/ExampleHeatmap2.png

.. tip:: **More examples** can be found in our `Gallery <http://deeptools.readthedocs.org/en/latest/content/example_gallery.html#normalized-chip-seq-signals-and-peak-regions>`_.

