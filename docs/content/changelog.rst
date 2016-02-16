Changes in deepTools2.0
========================

.. contents:: 
    :local:

Major changes
-------------

.. note:: The major changes encompass features for **increased efficiency**, **new sequencing data types**, and **additional plots**, particularly for QC.

Moreover, deepTools modules can now be used by other python programs.
The :ref:`api` is part of the new documentation.

Accommodating additional data types
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* correlation and comparisons can now be calculated for **bigWig files** (in addition to BAM files) using ``multiBigwigSummary`` and ``bigwigCompare``

* **RNA-seq:** split-reads are now natively supported
 
* **MNase-seq:** using the new option ``--MNase`` in ``bamCoverage``, one can now compute read coverage only taking the 2 central base pairs of each mapped fragment into account.

Structural updates
^^^^^^^^^^^^^^^^^^^

* All modules have comprehensive and automatic tests that evaluate proper functioning after any modification of the code.
* Virtualization for stability: we now provide a ``docker`` image and enable the easy deployment of deepTools via the Galaxy ``toolshed``.
* Our documentation is now version-aware thanks to readthedocs and ``sphinx``.
* The API is public and documented.

Renamed tools
^^^^^^^^^^^^^

* **heatmapper** to :doc:`tools/plotHeatmap`
* **profiler** to :doc:`tools/plotProfile`
* **bamCorrelate** to :doc:`tools/multiBamSummary`
* **bigwigCorrelate** to :doc:`tools/multiBigwigSummary`
* **bamFingerprint** to :doc:`tools/plotFingerprint`


Increased efficiency
^^^^^^^^^^^^^^^^^^^^

* We dramatically improved the **speed** of bigwig related tools (:doc:`tools/multiBigwigSummary` and ``computeMatrix``) by using the new `pyBigWig module <https://github.com/dpryan79/pyBigWig>`_.

* It is now possible to generate one composite heatmap and/or meta-gene image based on **multiple bigwig files** in one go (see :doc:`tools/computeMatrix`, :doc:`tools/plotHeatmap`, and :doc:`tools/plotProfile` for examples)

* ``computeMatrix`` now also accepts multiple input BED files. Each is treated as a group within a sample and is plotted independently.

* We added **additional filtering options for handling BAM files**, decreasing the need for prior filtering using tools other than deepTools: The ``--samFlagInclude`` and ``--samFlagExclude`` parameters can, for example, be used to only include (or exclude) forward reads in an analysis.

* We separated the generation of read count tables from the calculation of pairwise correlations that was previously handled by ``bamCorrelate``. Now, read counts are calculated first using ``multiBamSummary`` or ``multiBigWigCoverage`` and the resulting output file can be used for calculating and plotting pairwise correlations using ``plotCorrelation`` or for doing a principal component analysis using ``plotPCA``.

New features and tools
^^^^^^^^^^^^^^^^^^^^^^

* Correlation analyses are no longer limited to BAM files -- bigwig files are possible, too! (see :doc:`tools/multiBigwigSummary`)

* Correlation coefficients can now be computed even if the data contains NaNs.

* Added **new quality control** tools:
      - use :doc:`tools/plotCoverage` to plot the coverage over base pairs
      - use :doc:`tools/plotPCA` for principal component analysis
      - :doc:`tools/bamPEFragmentSize` can be used to calculate the average fragment size for paired-end read data
      
* Added the possibility for **hierarchical clustering**, besides *k*-means to ``plotProfile`` and ``plotHeatmap``

* ``plotProfile`` has many more options to make compelling summary plots


Minor changes
-------------

Changed parameters names and settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``computeMatrix`` can now read files with DOS newline characters.
* ``--missingDataAsZero`` was renamed to ``--skipNonCoveredRegions`` for clarity in ``bamCoverage`` and ``bamCompare``.
* Read extension was made optional and we removed the need to specify a default fragment length for most of the tools: ``--fragmentLength`` was thus replaced by the new optional parameter ``--extendReads``.
* Added option ``--skipChromosomes`` to ``multiBigwigSummary``, which can be used to, for example, skip all 'random' chromosomes.
* Added the option for adding titles to QC plots.

Bug fixes
^^^^^^^^^

* Resolved an error introduced by ``numpy version 1.10`` in ``computeMatrix``.
* Improved plotting features for ``plotProfile`` when using as plot type: 'overlapped_lines' and 'heatmap'
* Fixed problem with BED intervals in ``multiBigwigSummary`` and ``multiBamSummary`` that returned wrongly labeled raw counts.
* ``multiBigwigSummary`` now also considers chromosomes as identical when the names between samples differ by 'chr' prefix, e.g. chr1 vs. 1.
* Fixed problem with wrongly labeled proper read pairs in a BAM file. We now have additional checks to determine if a read pair is a proper pair: the reads must face each other and are not allowed to be farther apart than 4x the mean fragment length.
* For ``bamCoverage`` and ``bamCompare``, the behavior of ``scaleFactor`` was updated such that now, if given in combination with the normalization options (``--normalizeTo1x`` or ``--normalizeUsingRPKM``), the given scaling factor will be multiplied with the factor computed by the respective normalization method.


