# deepTools overview <a name="overview"></a>

deepTools consists of a set of modules that can be used independently to work with mapped reads. We have subdivided such tasks into *quality controls* (QC), *normalizations* and *visualizations*.

Here's a concise summary of the tools - if you would like more detailed information about the individual tools and example figures, follow the links in the table. 

| tool | type | input files | main output file(s) | application |
|------|------|-------------|---------------------|-------------|
| **bamCorrelate** | QC | 2 or more BAM | clustered heatmap | Pearson or Spearman correlation between read distributions |
| **bamFingerprint** | QC | 2 BAM | 1 diagnostic plot | assess enrichment strength of a ChIP sample |
| **computeGCbias** | QC | 1 BAM | 2 diagnostic plots | calculate the exp. and obs. GC distribution of reads|
| **correctGCbias** | QC | 1 BAM, output from computeGCbias | 1 GC-corrected BAM | obtain a BAM file with reads distributed according to the genome's GC content|
| **bamCoverage** | normalization | BAM | bedGraph or bigWig | obtain the normalized read coverage of a single BAM file |
| **bamCompare** | normalization | 2 BAM | bedGraph or bigWig | normalize 2 BAM files to each other using a mathematical operation of your choice (e.g. log2ratio, difference)|
| **computeMatrix** | visualization | 1 bigWig, 1 BED | gzipped table, to be used with heatmapper or profiler | compute the values needed for heatmaps and summary plots |
| **heatmapper** | visualization | computeMatrix output | heatmap of read coverages | visualize the read coverages for genomic regions |
| **profiler**| visualization | computeMatrix output | summary plot | visualize the average read coverages over a group of genomic regions | 