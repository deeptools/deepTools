`**WIKI-START** <Home>`__ > `**About deepTools** <About-deepTools>`__

-  `Why we built deepTools <#why>`__
-  `How we use deepTools <#howWe>`__
-  `What deepTools can do <#overview>`__

Why we built deepTools 
-----------------------

The main reason why deepTools was started is the simple fact that in
2011 we could not find tools that met all our needs for NGS data
analysis. While there were individual tools for separate tasks, we
wanted software that would fulfill *all* of the following criteria:

-  **efficiently extract reads from BAM files** and perform various
   computations on them
-  **turn BAM files of aligned reads into bigWig files** using different
   normalization strategies
-  make use of **multiple processors** (speed!)
-  generation of **highly customizable images** (change colours, size,
   labels, file format etc.)
-  enable **customized down-stream analyses** which requires that every
   data set that is being produced can be stored by the user
-  **modular approach** - compatibility, flexibility, scalability (i.e.
   we can add more and more modules making use of established methods)

| The flow chart below depicts the different tool modules that are
currently available within deepTools (deepTools modules are written in
bold red and black font). For more information on a typical analysis
pipeline, read the text below and `What deepTools can do <#overview>`__.
| |flowChartI|. If you the file names in the figure mean nothing to you,
please make sure to check our
`Glossary <https://github.com/fidelram/deepTools/wiki/Glossary>`__.

How we use deepTools 
---------------------

You will find many examples from ChIP-seq analyses in this tutorial, but
this does not mean that deepTools is restricted to ChIP-seq data
analysis. However, some tools, such as *bamFingerprint* specifically
address ChIP-seq-issues.

`Here <https://docs.google.com/file/d/0B8DPnFM4SLr2UjdYNkQ0dElEMm8/edit?usp=sharing>`__
are slides that we used for teaching at the University of Freiburg that
contain more details on the deepTools usage and aims.

| As shown in the flow chart above, our work usually begins with one or
more
`FASTQ <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-fastq>`__
file(s) of deeply-sequenced samples. After a first quality control using
`FASTQC <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`__,
we align the reads to the reference genome, e.g. using
`bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`__.
| We then use deepTools to assess the quality of the aligned reads:

#. **Correlation between
   `BAM <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bam>`__
   files** (*bamCorrelate*). This is a very basic test to see whether
   the sequenced and aligned reads meet your expectations. We use this
   check to assess the reproducibility - either between replicates
   and/or between different experiments that might have used the same
   antibody or the same cell type etc. For instance, replicates should
   correlate better than differently treated samples.
#. **GC bias check** (*computeGCbias*). Many sequencing protocols
   require several rounds of PCR-based amplification of the DNA to be
   sequenced. Unfortunately, most DNA polymerases used for PCR introduce
   significant GC biases as they prefer to amplify GC-rich templates.
   Depending on the sample (preparation), the GC bias can vary
   significantly and we routinely check its extent. In case we need to
   compare files with different GC biases, we use the *correctGCbias*
   module to match the GC bias.
   See the paper by `Benjamini and
   Speed <http://nar.oxfordjournals.org/content/40/10/e72>`__ for many
   insights into this problem.
#. **Assessing the ChIP strength**. We do this quality control to get a
   feeling for the signal-to-noise ratio in samples from ChIP-seq
   experiments. It is based on the insights published by `Diaz et
   al. <http://www.degruyter.com/view/j/sagmb.2012.11.issue-3/1544-6115.1750/1544-6115.1750.xml>`__.

Once we're satisfied by the basic quality checks, we normally **convert
the large
`BAM <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bam>`__
files into a leaner data format, typically
`bigWig <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bigwig>`__**.
bigWig files have several advantages over BAM files that mainly stem
from their significantly decreased size:

-  useful for data sharing & storage
-  intuitive visualization in Genome Browsers (e.g.
   `IGV <http://www.broadinstitute.org/igv/>`__)
-  more efficient downstream analyses are possible

The deepTools modules *bamCompare* and *bamCoverage* do not only allow
the simple conversion from BAM to bigWig (or
`bedGraph <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bedgraph>`__
for that matter), **the main reason why we developed those tools was
that we wanted to be able to *normalize* the read coverages** so that we
could compare different samples despite differences in sequencing depth,
GC biases and so on.

Finally, once all the files have passed our visual inspections, the fun
of downstream analyses with *heatmapper* and *profiler* can begin!

deepTools overview 
-------------------

deepTools consists of a set of modules that can be used independently to
work with mapped reads. We have subdivided such tasks into *quality
controls* (QC), *normalizations* and *visualizations*.

Here's a concise summary of the tools - if you would like more detailed
information about the individual tools and example figures, follow the
links in the table.

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

--------------

.. |flowChartI| image:: https://raw.github.com/fidelram/deepTools/master/examples/flowChart_BAMtoBIGWIG.png
