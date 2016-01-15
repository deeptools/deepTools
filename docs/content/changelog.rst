Changes in deepTools2.0
========================

Major changes
-------------

The major changes encompass features for **increased efficiency**, 
**new sequencing data types**, and **additional plots**.

Moreover, deepTools modules can now be used by other python programs.
The :ref:`api` is now part of the documentation.

One of the most visible changes is certainly the move of the
documentation from the previous github-hosted wiki to http://deeptools.readthedocs.org.


Additional data types
^^^^^^^^^^^^^^^^^^^^^^

* **RNA-seq:** split-reads are now natively supported
 
* **MNase-seq:** using the new option ``--MNase`` in ``bamCoverage``, one can now compute read coverage only taking the 2 central base pairs of each mapped fragment into account.
 

Increased efficiency
^^^^^^^^^^^^^^^^^^^^^

* We dramatically improved the **speed** of bigwig related tools (``bigwigCorrelate`` and ``computeMatrix``) by using the new `pyBigWig module <https://github.com/dpryan79/pyBigWig>`_.

* It is now possible to generate one composite heatmap and/or meta-gene image based on **multiple bigwig files** in one go (see :doc:`tools/computeMatrix`, :doc:`tools/plotHeatmap`, and :doc:`tools/plotProfile` for examples)

* ``computeMatrix`` also now accepts multiple input BED files. Each is treated as a group within a sample and is plotted independently.

* We added **additional filtering options for handling BAM files**, decreasing the need for prior filtering using tools other than deepTools: The ``--samFlagInclude`` and ``--samFlagExclude`` parameters can, for example, be used to only include (or exclude) forward reads in an analysis.

* We separated the generation of read count tables from the calculation of pairwise correlations that was previously handled by ``bamCorrelate``. Now, read counts are calculated first using ``multiBamCoverage`` or ``multiBigWigCoverage`` and the resulting output file can be used for calculating and plotting pairwise correlations using ``plotCorrelation`` or for doing a principal component analysis using ``plotPCA``.

New features and tools
^^^^^^^^^^^^^^^^^^^^^^^

* Correlation analyses are no longer limited to BAM files -- bigwig files are possible, too! (see :doc:`tools/bigwigCorrelate`)
* Correlation coefficients can now be computed even if the data contains NaNs.
* Added **new quality control** tools:
      - use :doc:`tools/plotCoverage` to plot the coverage over base pairs
      - use :doc:`tools/plotPCA` for principal component analysis
* Added the possibility for **hierarchical clustering**, besides *k*-means to ``plotProfile`` and ``plotHeatmap``


Minor changes
-------------

Changed parameter names and settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Renamed **heatmapper** to ``plotHeatmap`` and **profiler** to ``plotProfile``
* ``computeMatrix`` can now read files with DOS newline characters.
* ``--missingDataAsZero`` was renamed to ``--skipNonCoveredRegions`` for clarity in ``bamCoverage`` and ``bamCompare``.
* Read extension was made optional and we removed the need to specify a default fragment length for most of the tools: ``--fragmentLength`` was thus replaced by the new optional parameter ``--extendReads``.
* Added option ``--skipChromosomes`` to ``bigwigCorrelate``, for example to skip all 'random' chromosomes.
* Added the option for adding titles to QC plots.

Bug fixes
^^^^^^^^^^
* ``bigwigCorrelate`` now also considers chromosomes as identical when the names between samples differ by 'chr' prefix, e.g. chr1 vs. 1.
* Resolved an error introduced by numpy version 1.10 in ``computeMatrix``.
* Improved plotting features for ``tools/plotProfile`` when using as plot type: 'overlapped_lines' and 'heatmap'
* Fixed problem with BED intervals in ``bigwigCorrelate`` and ``bamCorrelate`` that returned wrongly labeled raw counts.
* Fixed problem with wrongly labeled proper read pairs in a BAM file. We now have additional checks to determine if a read pair is a proper pair.
* For ``bamCoverage`` and ``bamCompare``, behaviour of ``scaleFactor`` was updated such that now, if given in combination with the normalization options (``--normalizeTo1x`` or ``--normalizeUsingRPKM``), the given scaling factor will multiply the scale factor computed for the normalization methods.