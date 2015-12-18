Changes deepTools2.0
====================

Major changes
^^^^^^^^^^^^^

 * :doc:`tools/computeMatrix` now accepts multiple bigwig files that can later be plotted together as heatmaps
   one after the other or as multiple lines in the same plot. See the documentation of :doc:`tools/plotHeatmap`
   and :doc:`tools/plotProfile` for examples.

 * :doc:`tools/computeMatrix` also now accepts multiple input BED files. Each is treated as a group within a sample
   and is plotted independently.

 * Dramatically improved the speed of bigwig related tools (:doc:`tools/bigwigCorrelate` and :doc:`tools/computeMatrix`)
   by using the new `pyBigWig module <https://github.com/dpryan79/pyBigWig>`_.

 * Added support for split reads (most commonly found in RNA-seq data).

 * Added ``--samFlagInclude`` and ``--samFlagExclude`` parameters. This is useful to for example
   only include forward reads (or only reverse reads) in an analysis.

 * The documentation was migrated to http://deeptools.readthedocs.org

 * deepTools modules can now be used by other python programs. The :ref:`api` is now part of the documentation.

 * In this new release, most of the core code was rewriting to facilitate API usage and for optimization.

Minor changes
^^^^^^^^^^^^^

 * Added new quality control tool :doc:`tools/plotCoverage` to plot the coverage over base pairs for multiple samples
 * ``--missingDataAsZero`` was renamed to ``--skipNonCoveredRegions`` for clarity in :doc:`tools/bamCoverage`
   and :doc:`tools/bamCompare`.
 * Added new analysis tool :doc:`tools/plotPCA` to visualize the results of :doc:`tools/bamCorrelate`
   or :doc:`tools/bigwigCorrelate` using principal component analysis.
 * Added new option ``--MNase`` in :doc:`tools/bamCoverage` that computes reads coverage only considering two
   base pairs at the center of the fragment.
 * Read extension was made optional and removed the need to specify a default fragment length for most of the tools.
   and ``--fragmentLentgh parameters`` were replaced by the new optional parameter ``--extendReads``.
 * Renamed **heatmapper** to :doc:`tools/plotHeatmap` and **profiler** to :doc:`tools/plotProfile`
 * Added hierarchical clustering, besides *k*-means to :doc:`plotProfile` and :doc:`tools/plotHeatmap`
 * Improved plotting features for :doc:`tools/plotProfile` when using as plot type: 'overlapped_lines' and 'heatmap'
 * Resolved an error introduced by numpy version 1.10 in :doc:`tools/computeMatrix:
 * Plotting of correlations (from :doc:`tools/bamCorrelate` or :doc:`tools/bigwigCorrelate`) is now
   separated from the computation of the underlying data. A new tool, :doc:`tools/plotCorrelation` was added. This tool
   can plot correlations as heatmaps or as scatter plots and includes options to adjust a large array of visual features.
 * Fixed problem with bed intervals in :doc:`tools/bigwigCorrelate` and :doc:`tools/bamCorrelate` and a
   user specified region that returned wrongly labeled raw counts.
 * Correlation coefficients can now be computed even if the data contains NaNs.
 * :doc:`tools/computeMatrix` can now read files with DOS newline characters.
 * Added option ``--skipChromosomes`` to  :doc:`tools/bigwigCorrelate`, for example to skip all
   'random' chromosomes. :doc:`tools/bigwigCorrelate` now also considers chromosomes as identical
   when the names between samples differ by 'chr' prefix 'chr'. E.g. chr1 vs. 1
 * For :doc:`tools/bamCoverage` and :doc:`tools/bamCompare`, behaviour of scaleFactor was updated such that now,
   if given in combination with the normalization options (normalize to 1x or normalize using RPKM) the given scaleFactor
   will multiply the scale factor computed for the normalization methods.
 * Fixed problem with wrongly labeled proper pairs in a bam file. deepTools adds further checks to
   determine if a read pair is a proper pair.
 * Added titles to QC plots,
