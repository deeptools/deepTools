`How can I do...? Some exemplary protocols <#HowTo>`__
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  `I have downloaded/received a BAM file - how do I generate a file I
   can look at in a Genome Browser? <#FASTQ2IGV>`__

-  `How can I assess the reproducibility of my sequencing
   replicates? <#repCorr>`__

-  `How do I know whether my sample is GC biased? And if yes, how do I
   correct for it? <#GC>`__

-  `How do I get an input-normalized ChIP-seq coverage
   file? <#InputNorm>`__

-  `How can I compare the ChIP strength for different ChIP
   experiments? <#fprint>`__

-  `How do I get a (clustered) heatmap of sequencing-depth-normalized
   read coverages around the transcription start site of all
   genes? <#HM>`__

-  `How can I compare the average signal for X- and autosomal genes for
   2 or more different sequencing experiments <#multiprofiler>`__

`Go to the deepTools Galaxy help page <Galaxy>`__
                                                 

--------------

How can I do...? 
-----------------

This section is meant to give you quick guidance to specific tasks you
may want to perform. We're using screenshots from Galaxy here, if you're
using the command-line version, you can easily follow the given examples
by typing the program name and the help option (e.g.
/deepTools/bin/bamCoverage --help) which will show you all the
parameters and options, most of them named very similarly to those in
Galaxy.

For each "recipe" here, you will find the screenshot of the tool and the
input parameters on the left hand side (we marked non-default,
*user-specified entries*) and screenshots of the output on the right
hand side. Do let us know if you spot things that are missing, should be
explained better or are plain confusing!

There are many more ways in which you can use `deepTools
Galaxy <http://deeptools.ie-freiburg.mpg.de/>`__ than those described
here, so be creative once you're comfortable with using them. For
detailed explanations of what the tools do, follow the links.

    All recipes assume that you have uploaded your files into a Galaxy
    instance with a deepTools installation, e.g. `deepTools
    Galaxy <http://deeptools.ie-freiburg.mpg.de/>`__.

    If you would like to try out the protocols with **sample data**, go
    to `deepTools Galaxy <http://deeptools.ie-freiburg.mpg.de/>`__ →
    "Shared Data" → "Data Libraries" → "Sample Data". Use one file from
    each folder. E.g., import 1 BAM file from the folder "mapped reads",
    1 bigwig file from the folder "normalized read coverages"
    (preferably "Dmel\_log2ratio....bigwig", and the .bed file
    "Drosophila housekeeping genes" from the folder "annotation data"
    into your current Galaxy history. For testing our protocols via the
    command line, you can download the sample files to your computer by
    clicking on the triangle right next to the file name.

--------------

I have downloaded/received a `BAM <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bam>`__ file - how do I generate a file I can look at in a Genome Browser?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  tool:
   `bamCoverage <https://github.com/fidelram/deepTools/wiki/Normalizations>`__
-  input: your BAM file

Note: BAM files can also be viewed in Genome Browsers, however, they're
large and tend to freeze the applications. Generating bigWig files of
read coverages will help you a lot in this regard. In addition, if you
have more than one sample you'd like to look at, it is helpful to
normalize all of them to 1x sequencing depth.

| 
| 
| 

--------------

How can I assess the reproducibility of my sequencing replicates?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  tool:
   `bamCorrelate <https://github.com/fidelram/deepTools/wiki/QC>`__
-  input: BAM files

   -  you can compare as many samples as you want - the more you put at
      the same time, the longer the computation takes

-  output: heatmap of correlations - the closer two samples are to each
   other, the more similar their read coverages

| 
| 
| 

--------------

How do I know whether my sample is GC biased? And if yes, how do I correct for it?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  you need a BAM file of your sample in question
-  use the tool
   `computeGCbias <https://github.com/fidelram/deepTools/wiki/QC>`__ on
   that BAM file (default settings, just make sure your reference genome
   and genome size are matching)

| 
| 
| 

-  have a look at the image that is produced and compare it to the
   examples `here <https://github.com/fidelram/deepTools/wiki/QC>`__
-  if your sample shows an almost linear increase in exp/obs coverage
   (on the log scale of the lower plot), then you should consider
   correcting the GC bias - *if* you think that the biological
   interpretation of this data would otherwise be compromised (e.g. by
   comparing it to another sample that does not have an inherent GC
   bias)

   -  the GC bias can be corrected with the tool
      `correctGCbias <https://github.com/fidelram/deepTools/wiki/Normalizations>`__
      using the second output of the computeGCbias tool that you had to
      run anyway

   -  CAUTION!! correctGCbias will add reads to otherwise depleted
      regions (typically GC-poor regions), that means that you should
      **not** remove duplicates in any downstream analyses based on the
      GC-corrected BAM file (we therefore recommend to remove duplicates
      before doing the correction so that only those duplicate reads are
      kept that were produced by the GC correction procedure)

| 
| 
| 

--------------

How do I get an input-normalized ChIP-seq coverage file?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  input: you need two BAM files, one for the input, one for the
   ChIP-seq experiment
-  tool:
   `bamCompare <https://github.com/fidelram/deepTools/wiki/Normalizations>`__
   with ChIP = treatment, input = control sample

| 
| 
| 

--------------

How can I compare the ChIP strength for different ChIP experiments?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  tool:
   `bamFingerprint <https://github.com/fidelram/deepTools/wiki/QC>`__
-  input: as many BAM files as you'd like to compare. Make sure you get
   all the labels right!

| 
| 
| 

--------------

How do I get a (clustered) heatmap of sequencing-depth-normalized read coverages around the transcription start site of all genes?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  tools:
   `computeMatrix <https://github.com/fidelram/deepTools/wiki/Visualizations>`__,
   then
   `heatmapper <https://github.com/fidelram/deepTools/wiki/Visualizations>`__
-  inputs:

   -  1 bigWig file of normalized read coverages (e.g. the result of
      bamCoverage or bamCompare)
   -  1 BED or INTERVAL file of genes, e.g. obtained through Galaxy via
      "Get Data" → "UCSC main table browser" → group: "Genes and Gene
      Predictions" → (e.g.) "RefSeqGenes" → send to Galaxy (see
      screenshots below)

| 
| 
| 

-  use
   `computeMatrix <https://github.com/fidelram/deepTools/wiki/Visualizations>`__
   with the bigWig file and the BED file
-  indicate "reference-point" (and whatever other option you would like
   to tune, see screenshot below)

| 
| 
| 

-  use the output from computeMatrix with
   `heatmapper <https://github.com/fidelram/deepTools/wiki/Visualizations>`__

   -  if you would like to cluster the signals, choose "kmeans
      clustering" (last option of "advanced options") with a reasonable
      number of clusters (usually between 2 to 7)

| 
| 
| 

--------------

How can I compare the average signal for X- and autosomal genes for 2 or more different sequencing experiments?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure you're familiar with computeMatrix and profiler before using
this protocol.

-  tools:

   -  Filter data on any column using simple expressions
   -  computeMatrix
   -  profiler
   -  (plotting the summary plots for multiple samples)

-  inputs:

   -  several bigWig files (one for each sequencing experiment you would
      like to compare)
   -  two BED files, one with X-chromosomal and one with autosomal genes
      \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_

How to obtain a BED file for X chromosomal and autosomal genes each
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

#. download a full list of genes via "Get Data" → "UCSC main table
   browser" → group:"Genes and Gene Predictions" → tracks: (e.g.)
   "RefSeqGenes" → send to Galaxy

#. filter the list twice using the tool **"Filter data on any column
   using simple expressions"**

   -  first use the expression: c1=="chrX" to filter the list of all
      genes → this will generate a list of X-linked genes
   -  then re-run the filtering, now with c1!="chrX" which will generate
      a list of genes that do not belong to chromosome X (!= indicates
      "not matching")

Compute the average values for X and autosomal genes
''''''''''''''''''''''''''''''''''''''''''''''''''''

-  use
   `computeMatrix <https://github.com/fidelram/deepTools/wiki/Visualizations>`__
   for **each** signal file (bigWig) (you only need to specify all the
   parameters once, then use the re-run button underneath the first data
   set and just replace the signal file with the next one)

   -  supply both filtered BED files (click on "Add new regions to plot"
      once) and label them
   -  indicate the corresponding signal file
   -  make sure to **re-name** every data set in the history once
      computeMatrix is done so that you can easily keep track of which
      matrix was based on which bigWig file (you can always find these
      information by clicking on the i-button in the respective data
      set)

-  now use
   `profiler <https://github.com/fidelram/deepTools/wiki/Visualizations>`__
   for every file you generated with computeMatrix

   -  important: display the "advanced output options" and select "save
      the data underlying the average profile" → this will generate a
      table in addition to the summary plot images

-  now you have at least 2 separate images of profiles - one for each
   bigWig file - you can either leave it like this or use another script
   that will plot all the summary plots in one image at once

   -  this tool is called "Plotting the summary plots for multiple
      signals"
   -  it uses the tables generated by profiler
   -  for each group of genes (in this case, X and autosomal genes = 2
      groups), you can assign a color

The result could look like this:

| 
| 
| 

As you may have noticed, this task requires several steps that are
repeated. Here is a screenshot of how the Galaxy workflow would look
like (you can find it under "Shared Data" → "Published Workflows" →
"Summary plots for X and autosomal genes" where we have constructed it
with the example histone marks from the Data Library. Be aware that
running this workflow will take up quite some computation timing, but it
won't require much input from your part - so start if before you go off
for lunch...)

If you're not sure how to use the published workflow, please read `this
entry <#workflow>`__ or go to the central `Galaxy learning page full of
tutorials <https://wiki.galaxyproject.org/Learn/Screencasts#Tutorials>`__.

| 
| 
| 

--------------

