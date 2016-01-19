Which tools can I find in the deepTools Galaxy?
-----------------------------------------------

As mentioned before, each Galaxy installation can be tuned to the
individual interests.
Our goal is to provide a Galaxy that enables you to **quality check, process and normalize and subsequently visualize your data obtained by high-throughput DNA sequencing**.

.. tip:: If you do not know the difference between a BAM and a BED file, that's fine. You can read up on them in our :doc:`help_glossary`.

.. tip:: For more specific help, check our :doc:`help_faq_galaxy` and the :doc:`example_step_by_step`.

We provide the following kinds of tools:

.. contents:: 
    :local:

deepTools
^^^^^^^^^^

The most important category is called **"deepTools"** that contains all the main tools we have developed.

Tools for BAM and bigWig file processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+------------------------+--------------------------------------------------------------------------------+
| ``multiBamCoverage``   | get read counts for the binned genome or user-specified regions                |
+------------------------+--------------------------------------------------------------------------------+
| ``multiBigwigSummary`` | calculate score summaries for the binned genome or user-specified regions      |
+------------------------+--------------------------------------------------------------------------------+
| ``correctGCBias``      | obtain a BAM file with reads distributed according to the genome's GC content  |
+------------------------+--------------------------------------------------------------------------------+
| ``bamCoverage``        | obtain the normalized read coverage of a single BAM file                       |
+------------------------+--------------------------------------------------------------------------------+
| ``bamCompare``         | normalize 2 BAM files to each other (e.g. log2ratio, difference)               |
+------------------------+--------------------------------------------------------------------------------+
| ``bigwigCompare``      | normalize the scores of two bigWig files to each other (e.g., ratios)          |
+------------------------+--------------------------------------------------------------------------------+
| ``computeMatrix``      | compute the values needed for heatmaps and summary plots                       |
+------------------------+--------------------------------------------------------------------------------+

Tools for QC
^^^^^^^^^^^^^

+-----------------------+-------------------------------------------------------------------------------------------------------+
| ``plotCorrelation``   | calculate and visualize the pairwise Spearman or Pearson correlation of read counts (or other scores) |
+-----------------------+-------------------------------------------------------------------------------------------------------+
| ``plotPCA``           | perform PCA and visualize the results                                                                 |
+-----------------------+-------------------------------------------------------------------------------------------------------+
| ``plotFingerprint``   | assess the ChIP enrichment strength                                                                   |
+-----------------------+-------------------------------------------------------------------------------------------------------+
| ``bamPEFragmentSize`` | obtain the average fragment length for paired-end samples                                             |
+-----------------------+-------------------------------------------------------------------------------------------------------+
| ``computeGCBias``     | assess the GC bias by calculating the expected and observed GC distribution of aligned reads          |
+-----------------------+-------------------------------------------------------------------------------------------------------+
| ``plotCoverage``      | obtain the normalized read coverage of a single BAM file                                              |
+-----------------------+-------------------------------------------------------------------------------------------------------+

Heatmaps and summary plots
^^^^^^^^^^^^^^^^^^^^^^^^^^

+-------------------+-------------------------------------------------------------------------------------------+
| ``plotHeatmap``   | visualize read counts or other scores in heatmaps with one row per genomic region         |
+-------------------+-------------------------------------------------------------------------------------------+
| ``plotProfile``   | visualize read counts or other scores using average profiles (e.g., meta-gene profiles)   |
+-------------------+-------------------------------------------------------------------------------------------+


We have compiled several sources of detailed information specifically
about the usage of deepTools:

1. General overview of `how we use deep Tools <About-deepTools>`__
2. Each individual tool is described in more detail on separate pages -
   just follow the links in the table above
3. For each tool, you will find specific explanations within the
   `deepTools Galaxy <https://urldefense.proofpoint.com/v2/url?u=http-3A__deeptools.ie-2Dfreiburg.mpg.de_&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=1xpNM-7I4Z6-ZIJErjnO726mjBKdGH92RCWOc5kGh-U&e= >`__ main
   frame, too.
4. the `example workflows <Example-workflows>`__ might help to get a
   feeling for the kinds of analyses than can be done with `deepTools
   Galaxy <https://urldefense.proofpoint.com/v2/url?u=http-3A__deeptools.ie-2Dfreiburg.mpg.de_&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=1xpNM-7I4Z6-ZIJErjnO726mjBKdGH92RCWOc5kGh-U&e= >`__

Peak calling
^^^^^^^^^^^^^^

In ChIP-seq analysis, peak calling algorithms are essential downstream
analysis tools to identify regions of significant enrichments (i.e.
where the ChIP sample contained significantly more sequenced reads than
the input control sample). By now, there must be close to 100 programs
out there (see `Wilbanks et
al. <https://urldefense.proofpoint.com/v2/url?u=http-3A__www.plosone.org_article_info-253Adoi-252F10.1371-252Fjournal.pone.0011471&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=lhLQ7qst_E5ZweBT_PdS_mJIE9biseGu2DTBPk2papM&e= >`__
for a comparison of peak calling programs).

In contrast to deepTools that were developed for handling and generating
*continuous* genome-wide profiles, peak calling will result in a *list
of genomic regions*. Have a look at the screenshot to understand the
difference.

We have included the peak callers
`MACS <http://www.ncbi.nlm.nih.gov/pubmed/22936215>`__ and
`SICER <https://urldefense.proofpoint.com/v2/url?u=http-3A__bioinformatics.oxfordjournals.org_content_25_15_1952.full&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=4ZEmdk9_IT-qF0ZDdKMF6Z-vWNUrYB3r76ucGWLaCYo&e= >`__
within our Galaxy instance with
`MACS <http://www.ncbi.nlm.nih.gov/pubmed/22936215>`__ being the most
popular peak calling algorithm for the identification of localized
transcription factor binding sites while
`SICER <https://urldefense.proofpoint.com/v2/url?u=http-3A__bioinformatics.oxfordjournals.org_content_25_15_1952.full&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=4ZEmdk9_IT-qF0ZDdKMF6Z-vWNUrYB3r76ucGWLaCYo&e= >`__
was developed for diffuse ChIP-seq signals. Note that MACS version 1.14
is quite different from MACS version 2 (which has still not been
released officially).

Working with genomic intervals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Galaxy has 2 file formats to store lists of genomic regions:

-  INTERVAL

   -  tab-separated
   -  requirements:

      1. Column: chromosome
      2. Column: start position
      3. Column: end position

   -  all other columns can contain any value or character

-  BED

   -  very similar to INTERVAL, but stricter when it comes to what is
      expected to be kept in which column:

      -  

         1. to 3. Column: same as interval

      -  Column 4: name
      -  Column 5: score
      -  Column 6: strand

In case you would like to work with several lists of genomic regions,
e.g. generate a new list of regions that are found in two different
files etc., there are two categories of tools dedicated to performing
these tasks: \* Operate on genomic intervals \* BEDtools

Each tool's function is explained within Galaxy. Do browse those tools
as they will give you a very good glimpse of the scope of possible
analyses!

 #### Working with text files and tables In addition to deepTools that
were specifically developed for the handling of NGS data, we have
incorporated several standard Galaxy tools that enable you to manipulate
tab-separated files such as gene lists, peak lists, data matrices etc.

There are 3 main categories:

-  **Text manipulation**

   -  unlike Excel where you can easily interact with your text and
      tables via the mouse, data manipulations within Galaxy are
      strictly based on commands. If you feel like you would like to do
      something to certain *columns* of a data set, go through the tools
      of this category
   -  e.g. adding columns, cutting columns, pasting two files side by
      side, selecting random lines etc.
   -  a very useful tool of this category is called *Trim* - if you need
      to remove some characters from a column, this tool's for you! (for
      example, sometimes you need to adjust the chromosome naming
      between two files from different source - using *Trim*, you can
      remove the "chr" infront of the chromosome name)

-  **Filter and Sort**

   -  in addition to the common sorting and filtering, there's the very
      useful tool to *select lines that match an expression* (for
      example, using the expression *c1=='chrM'* will select all rows
      from a BED file with regions located on the mitochondrial
      chromosome)

-  **Join, Subtract, Group**
-  this category is very useful if you have several data sets that you
   would like to work with, e.g. by comparing them

**More help**

.. hint:: If you encounter a failing data set (marked in red), please send a bug report via the Galaxy bug report button and we will get in touch if you indicate your email address.

+-------------------------------------------------------------------------------+-----------------------------------------------------------------+
| `http://wiki.galaxyproject.org/Learn <http://wiki.galaxyproject.org/Learn>`_  | Help for Galaxy usage in general                                |
+-------------------------------------------------------------------------------+-----------------------------------------------------------------+
| `deepTools Galaxy FAQs <Galaxy-related-FAQs>`_                                | Frequently encountered issues with our specific Galaxy instance |
+-------------------------------------------------------------------------------+-----------------------------------------------------------------+
| deeptools@googlegroups.com                                                    | For issues not addressed in the FAQs                            |
+-------------------------------------------------------------------------------+-----------------------------------------------------------------+