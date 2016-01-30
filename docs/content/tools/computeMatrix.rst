computeMatrix
=============

.. contents:: 
    :local:

.. argparse::
   :ref: deeptools.computeMatrix.parse_arguments
   :prog: computeMatrix

Details
^^^^^^^^^^^^^^^

``computeMatrix`` is tightly connected to ``plotHeatmap`` and ``plotProfile``: it takes the values of all the signal files and all genomic regions that you would like to plot and computes the corresponding data matrix.

See :doc:`plotHeatmap` and :doc:`plotProfile` for example plots.

.. image:: ../../images/computeMatrix_overview.png

``computeMatrix`` has two main modes of use:

* for computing the signal distribution **relative to a point** (``reference-point``), e.g., the beginning or end of each genomic region
* for computing the signal **over a set of regions** (``scale-regions``) where all regions are scaled to the same size

.. tip:: ``computeMatrix`` can use multiple threads (``-p`` option), which significantly decreases the time for calculating the values.

In addition to generating the intermediate, gzipped file for ``plotHeatmap`` and ``plotProfile``, ``computeMatrix`` can also be used to simply output the values underlying the heatmap or to **filter and sort BED files** using, for example, the ``--skipZeros`` and the ``--sortUsing`` parameters.

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

