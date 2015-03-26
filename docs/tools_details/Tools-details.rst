`**WIKI-START** <Home>`__ > `**Tools overview** <Tools-details>`__

deepTools consists of a set of modules that can be used independently to
work with mapped reads. We have subdivided such tasks into *quality
controls*, *normalizations* and *visualizations*.

This table gives an overview of the tools that are available within the
current deepTools release. Most likely, we will add more modules in the
future.

For more detailed information, follow the links in the table or these
ones:

-  `tools for quality controls of aligned reads <QC>`__
-  `tools for bigWig generation <Normalizations>`__
-  `tools for visualization <Visualizations>`__

+---------------------------------------------------------+-----------------+------------------------------------+-------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+
| tool                                                    | type            | input files                        | main output file(s)                                   | application                                                                                                      |
+=========================================================+=================+====================================+=======================================================+==================================================================================================================+
| `bamCorrelate <QC#wiki-bamCorrelate>`__                 | QC              | 2 or more BAM                      | clustered heatmap                                     | Pearson or Spearman correlation between read distributions                                                       |
+---------------------------------------------------------+-----------------+------------------------------------+-------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+
| `bamFingerprint <QC#wiki-bamFingerprint>`__             | QC              | 2 BAM                              | 1 diagnostic plot                                     | assess enrichment strength of a ChIP sample                                                                      |
+---------------------------------------------------------+-----------------+------------------------------------+-------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+
| `computeGCbias <QC#wiki-computeGCbias>`__               | QC              | 1 BAM                              | 2 diagnostic plots                                    | calculate the exp. and obs. GC distribution of reads                                                             |
+---------------------------------------------------------+-----------------+------------------------------------+-------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+
| `correctGCbias <QC#wiki-correctGCbias>`__               | QC              | 1 BAM, output from computeGCbias   | 1 GC-corrected BAM                                    | obtain a BAM file with reads distributed according to the genome's GC content                                    |
+---------------------------------------------------------+-----------------+------------------------------------+-------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+
| `bamCoverage <Normalizations#wiki-bamCoverage>`__       | normalization   | BAM                                | bedGraph or bigWig                                    | obtain the normalized read coverage of a single BAM file                                                         |
+---------------------------------------------------------+-----------------+------------------------------------+-------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+
| `bamCompare <Normalizations#wiki-bamCompare>`__         | normalization   | 2 BAM                              | bedGraph or bigWig                                    | normalize 2 BAM files to each other using a mathematical operation of your choice (e.g. log2ratio, difference)   |
+---------------------------------------------------------+-----------------+------------------------------------+-------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+
| `computeMatrix <Visualizations#wiki-computeMatrix>`__   | visualization   | 1 bigWig, 1 BED                    | zipped file, to be used with heatmapper or profiler   | compute the values needed for heatmaps and summary plots                                                         |
+---------------------------------------------------------+-----------------+------------------------------------+-------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+
| `heatmapper <Visualizations#wiki-heatmapper>`__         | visualization   | computeMatrix output               | heatmap of read coverages                             | visualize the read coverages for genomic regions                                                                 |
+---------------------------------------------------------+-----------------+------------------------------------+-------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+
| `profiler <Visualizations#wiki-profiler>`__             | visualization   | computeMatrix output               | summary plot ("meta-profile")                         | visualize the average read coverages over a group of genomic regions                                             |
+---------------------------------------------------------+-----------------+------------------------------------+-------------------------------------------------------+------------------------------------------------------------------------------------------------------------------+

