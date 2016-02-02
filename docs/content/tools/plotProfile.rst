plotProfile
===========

.. argparse::
   :ref: deeptools.plotProfile.parse_arguments
   :prog: plotProfile

Details
^^^^^^^^

Like :doc:`plotHeatmap`, ``plotProfile`` simply takes the compressed matrix produced by ``computeMatrix`` and turns it into summary plots.

In addition to a large range of parameters for optimizing the visualization, you can also export the values underlying the profiles as tables.

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

Usage example
^^^^^^^^^^^^^^

The following example plots the signal profile over hg19 transcripts for our test ENCODE datasets. Note that the matrix contains multiple groups of regions (in this case, one for each present chromosome).

.. code:: bash

   $ plotProfile -m matrix_two_groups.gz \
        -out ExampleProfile1.png \
        --plotTitle "Test data profile"

.. image:: ../../images/test_plots/ExampleProfile1.png
<<<<<<< HEAD

``plotProfile`` has many options, including the ability to change the type of lines plotted and to plot by group rather than sample.
=======
>>>>>>> master

Here's the same data set, but plotted with a different set of parameters.

.. code:: bash

   $ plotProfile -m matrix_two_groups.gz \
        -out ExampleProfile2.png \
        --plotType=fill \ # add color between the x axis and the lines
        --perGroup \ # make one image per BED file instead of per bigWig file
        --colors red orange yellow green blue purple \
        --plotTitle "Test data profile"

.. image:: ../../images/test_plots/ExampleProfile2.png


