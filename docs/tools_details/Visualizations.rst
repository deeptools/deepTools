The modules for visualizing scores contained in
`bigWig <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bigwig>`__
files are separated into 1 tool that calculates the values
(*computeMatrix*) and 2 tools that contain many, many options to
fine-tune the plots (*heatmapper* and *profiler*). In other words:
computeMatrix generates the values that are the basis for heatmapper and
profiler.

| 
| 
| 

computeMatrix
-------------

This tool summarizes and prepares an intermediary file containing scores
associated with genomic regions that can be used afterwards to plot a
heatmap or a profile.

Genomic regions can really be anything - genes, parts of genes, ChIP-seq
peaks, favorite genome regions... as long as you provide a proper file
in
`BED <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bed>`__
or INTERVAL format. This tool can also be used to filter and sort
regions according to their score.

As indicated in the plot above, computeMatrix can be run with either one
of the two modes: **scaled regions** or **reference point**.

Please see the `example figures <#examples>`__ down below for
explanations of parameters and options.

Output files
^^^^^^^^^^^^

-  **obligatory**: zipped matrix of values to be used with *heatmapper*
   and/or *profiler*
-  **optional** (can also be generated with heatmapper or profiler in
   case you forgot to produce them in the beginning):

   -  BED-file of the regions sorted according to the calculated values
   -  list of average values per genomic bin
   -  matrix of values per genomic
      `bin <https://github.com/fidelram/deepTools/wiki/Glossary#terminology>`__
      per genomic interval

heatmapper
----------

| The heatmapper depicts values extracted from the
`bigWig <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bigwig>`__
file for each genomic region individually.
| It requires the output from computeMatrix and most of its options are
related to tweaking the visualization only. The values calculated by
computeMatrix are not changed.

Definitely check the examples at the bottom of the page to get a feeling
for how many things you can tune. In our
`Gallery <https://github.com/fidelram/deepTools/wiki/Gallery>`__ of
images generated with deepTools, you might find additional insights into
the possible use-cases of heatmapper and profiler.

profiler
--------

This tool plots the average enrichments over all genomic regions
supplied to computeMarix. It is a very useful complement to the
heatmapper, especially in cases when you want to compare the scores for
many different groups. Like heatmapper, profiler does not change the
values that were compute by computeMatrix, but you can choose between
many different ways to color and display the plots.

Example figures
---------------

Here you see a typical, not too pretty example of a heatmap. We will use
this example to explain several features of computeMatrix and
heatmapper, so do take a closer look.

1st example: Heatmap with all genes scaled to the one size and user-specified groups of genes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| 
| 
| 

As you can see, all genes have been scaled to the same size and the
(mean) values per
`bin <https://github.com/fidelram/deepTools/wiki/Glossary#terminology>`__
size (10 bp) are colored accordingly. In addition to the gene bodies, we
added 500 bp up- and down-stream of the genes.

The plot was produced with the following commands:

::

    $ /deepTools-1.5.2/bin/computeMatrix scale-regions --regionsFileName Dm.genes.indChromLabeled.bed \
    --scoreFileName PolII.bw --beforeRegionStartLength 500 --afterRegionStartLength 500 \
    --regionBodyLength 1500 --binSize 10 \
    --outFileName PolII_matrix_scaledGenes --sortRegions no

    $ /deepTools-1.5.2/bin/heatmapper --matrixFile PolII_matrix_scaledGenes \
    --outFileName PolII_indChr_scaledGenes.pdf \
    --plotTitle "Pol II" --whatToShow "heatmap and colorbar"

This is what you would have to select to achieve the same result within
Galaxy (pay attention to the fact that you will have to use two tools,
computeMatrix and heatmapper):

computeMatrix
             

| 
| 
| 

heatmapper
          

| 
| 
| 

The main difference between computeMatrix usage on the command line and Galaxy: the input of the regions file (BED)
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Note that we supplied just *one* BED-file via the command line whereas
in Galaxy we indicated three different files (one per chromosome).

| On the command line, the program expects a BED file where different
groups of genomic regions are concatenated into one file, where each
group should be indicated by "#group name" *following the last region*
of a particular group.
| The BED-file that was used here, contained 3 such lines and could be
prepared as follows:

::

     $ grep ^chr2 AllGenes.bed > Dm.genes.indChromLabeled.bed
     $ echo "#chr2" >> Dm.genes.indChromLabeled.bed
     $ grep ^chr3 AllGenes.bed >> Dm.genes.indChromLabeled.bed
     $ echo "#chr3" >> Dm.genes.indChromLabeled.bed
     $ grep ^chrX AllGenes.bed >> Dm.genes.indChromLabeled.bed
     $ echo "#chrX" >> Dm.genes.indChromLabeled.bed
     

In Galaxy, you can simply generate three different data sets starting
from a whole genome list of *Drosophila melanogaster* genes by using the
"Filter" tool ("Filter and Sort" --> "Filter") on the entries in the
first column three times:

#. c1=="chr2" --> Dm.genes.chr2.bed
#. c1=="chr3" --> Dm.genes.chr3.bed
#. c1=="chrX" --> Dm.genes.chrX.bed

Important parameters for optimizing the visualization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. **sorting of the regions**: The default of heatmapper is to sort the
   values in descending order. You can change that to ascending, no
   sorting at all or according to the size of the region (Using the
   ``--sort`` option on the command line or advanced options in Galaxy).
   We strongly recommend to leave the sorting option at "no sorting" for
   the initial computeMatrix step.
#. **coloring**: The default coloring by heatmapper is done using the
   python color map "RdYlBu", but this can be changed (--colorMap on the
   command line, advanced options within Galaxy).
#. **dealing with missing data**: You have certainly noticed that some
   gene bodies are depicted as white lines within the otherwise colorful
   mass of genes. Those regions are due to genes that, for whatever
   reason, did not have any read coverage in the bigWig file. There are
   several ways to handle these cases:

   -  **--skipZeros** this is useful when your data actually has a quite
      nice coverage, but there are 2 or 3 regions where you deliberately
      filtered out reads or you don't expect any coverage (e.g. hardly
      mapable regions). This will only work if the entire region does
      not contain a single value.
   -  **--missingDataAsZero** this option allows computeMatrix do
      interpret missing data points as zeroes. Be aware of the changes
      to the average values that this might cause.
   -  **--missingDataColor** this is in case you have very sparse data
      or were missing values make sense (e.g. when plotting methylated
      CpGs - half the genome should have no value). This option then
      allows you to pick out your favorite color for those regions. The
      default is black (was white when the above shown image was
      produced).

2nd example: Summary plots with all genes scaled to the one size and user-specified groups of genes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here's the **profiler** plot corresponding to the heatmap above. There's
one major difference though - do you spot it?

| 
| 
| 

We used the same
`BED <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bed>`__
file(s) as for the heatmap, hence the 3 different groups (1 per
chromosome). However, this time we used computeMatrix not with
*scale-regions* but with *reference-point* mode.

::

    $ /deepTools-1.5.2/bin/computeMatrix reference-point --referencePoint TSS \
    --regionsFileName Dm.genes.indChromLabeled.bed --scoreFileName PolII.bw \
    --beforeRegionStartLength 1000 --afterRegionStartLength 1000 \
    --binSize 10 --outFileName PolII_matrix_indChr_refPoint \
    --missingDataAsZero --sortRegions no

    $ /deepTools-1.5.2/bin/profiler --matrixFile PolII_matrix_indChr_refPoint \
    --outFileName profile_PolII_indChr_refPoint.pdf
    --plotType fill --startLabel "TSS" \
    --plotTitle "Pol II around TSS" --yAxisLabel "mean Pol II coverage" \
    --onePlotPerGroup

When you compare the profiler commands with the heatmapper commands, you
also notice that we made use of many more labeling options here, e.g.
``--yAxisLabel`` and a more specific title via ``-T``

This is how you would have obtained this plot in Galaxy (only the part
that's *different* from the above shown command for the scale-regions
version is shown):

computeMatrix
'''''''''''''

| 
| 
| 

profiler
''''''''

| 
| 
| 

3rd example: Heatmap with all genes scaled to the one size and kmeans clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Instead of supplying groups of regions on your own, you can use the
clustering function of heatmapper to get a first impression whether the
signal of your experiment can be easily clustered into two or more
groups of similar signal distribution.

Have a look at this example with two clusters. The values correspond to
log2ratios(ChIP/input) from a ChIP-seq experiment for RNA Polymerase II
in *Drosophila melanogaster*:

| 
| 
| 

The plot was produced with the following commands:

::

    $ /deepTools-1.5.2/bin/computeMatrix reference-point \
    --regionsFilenName Dm.genes.indChromLabeled.bed \
    --scoreFileName PolII.bw \
    --beforeRegionStartLength 500 --afterRegionStartLength 5000 \
    --binSize 50 \
    --outFileName PolII_matrix_TSS

    $ /deepTools-1.5.2/bin/heatmapper --matrixFile PolII_matrix_TSS \
    --kmeans 2 \
    --sortUsing region_length \
    --outFileName PolII_two_clusters.pdf \
    --plotTitle "Pol II" --whatToShow "heatmap and colorbar"

In Galaxy, these are the screenshots from the commands for computeMatrix
and heatmapper:

computeMatrix
             

| 
| 
| 

heatmapper
          

| 
| 
| 

| When the ``--kmeans`` option is chosen and more than 0 clusters are
specified, heatmapper will run the
`k-means <http://en.wikipedia.org/wiki/K-means_clustering>`__ clustering
algorithm. In this example, we wanted to divide *Drosophila
melanogaster* genes into two clusters. As you can see above, the
algorithm nicely identified two groups - one with mostly those genes
with lots of Pol II at the promoter region (top) from those genes
without Poll II at the promoter (bottom).
| Please note that the clustering will only work if the initial BED-file
used with computeMatrix contained only *one* group of genes.

The genes belonging to each cluster can be obtained by via
``--outFileSortedRegions`` on the command line and "advanced output
options in Galaxy". On the command line, this will result in a BED file
where the groups are separated by a hash tag. In Galaxy, you will obtain
individual data sets per cluster.

To have a better control on the clustering it is recommended to load the
matrix raw data into **specialized software like
`cluster3 <http://bonsai.hgc.jp/~mdehoon/software/cluster/>`__ or
`R <http://www.r-project.org/>`__**. You can obtain the matrix via the
option ``--outFileNameMatrix`` on the command line and by the "advanced
output options" in Galaxy. The order of the rows is the same as in the
output of the ``--outFileSortedRegions`` BED file.

--------------

