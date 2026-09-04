Changes in deepTools4.0
=======================

Plotting
--------

* Plots in general:
    - Removed Plotly for all graphics.
	- Optional ggplot theme (--ggplot).
	- Using standard Matplotlib color scheme.

* plotHeatmap and plotProfile:
	- Adjusted the legend position to be outside the plot area, preventing overlap with the subplots.
	- Split and wrapped longer sample names to avoid overlap with other labels.
	- Eliminated overlapping x and y-axis ticks for improved readability.
	- Removed the file extensions from input files when using them as labels in the plots.

* plotEnrichment:
	- Added spacing between bars in bar plots for improved visual clarity.
	- Retained only sample names, when provided with the complete file path and name.

* plotPCA:
    - Re-implemented with a lightweight scipy/numpy SVD backend (mirroring scikit-learn's PCA), instead of the previous plotly-based implementation.
	- New option to add labels for each point (--addLabels).
	- Expander for colors and markers, for example ``--colors red:3 blue:3`` is expanded as ``[red, red, red, blue, blue, blue]``.
	- Scree plot is showing lines for individual and accumulated variation.
	- Points are by default rainbow colored circles.

Core
----

* bamCoverage
    - --no_collapse flag to not merge bins with equal coverage values together.
 
* computeMatrix
    - Sorting ascend / descend no longer has subsorting by position.
    - --quiet / -q option no longer exists.
    - bed files in computeMatrix no longer support '#' to define groups.
    - 'chromosome matching' i.e. chr1 <-> 1, chrMT <-> MT is no longer performed.
	- metagene mode erroneously 'nan'ed the before and after values (if they fell outside of the feature). This is fixed now.
	- Rounding bahvior in matrix output only two decimals now, unscaled 5 and unscaled 3 prime are now strictly separated from the rest of the scaled region (for value calculation).
	
* normalization
    - Exactscaling is no longer an option, it's always performed.
	- SES option in bamCompare mode is no longer available.
	- blackList filtering is now performed on a position-based level. Meaning reads that overlap partially with the blacklist can still contribute to the signal.

* alignmentSieve
    - --genomeChunkLength is removed; --label, --smartLabels and --ignoreDuplicates are still available.
    - output order now exactly matches input order.

* new Rust-backed core
    - bamCoverage, bamCompare, computeMatrix, alignmentSieve and multiBamSummary are now backed by a Rust core (rayon multithreading, bigtools-based bigWig I/O), replacing the pure-Python implementations.
    - the previous Python implementations remain available as ``bamCoverage_old``, ``bamCompare_old``, ``computeMatrix_old``, ``alignmentSieve_old`` and ``multiBamSummary_old`` during the transition, but will be removed in a future release.
    - gzipped GTF/BED region files are now supported in multiBamSummary, computeMatrix and alignmentSieve; gzipped blacklist files are supported in all five Rust-backed tools.
    - large scale values have slightly different precision than before (f32 vs f64 internally).

Testing
-------
