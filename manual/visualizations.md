Visualization
==============

The modules for visualizing scores contained in [bigWig][] files are separated into a tool that calculates the values
(_computeMatrix_) and two tools that contain many, many options to fine-tune the plot (_heatmapper_ and _profiler_).

![flowChartII](https://raw.github.com/fidelram/deepTools/master/examples/flowChart_computeMatrixetc.png "Relationship between computeMatrix, heatmapper and profiler")

## Table of content

  * [computeMatrix](#computeMatrix)
  * [heatmapper](#heatmapper)
  * [profiler](#profiler)
  * [Example figures](#examples)

<a name="computeMatrix"/>
## computeMatrix
This tool summarizes and prepares an intermediary file containing
scores associated with genomic regions that can be used afterwards to
plot a heatmap or a profile. 

Genomic regions can really be anything - genes, parts of genes, ChIP-seq peaks, favorite genome regions...
as long as you provide a proper file in [BED][] or INTERVAL format. This tool can also be used to filter and sort regions
according to their score.

As indicated in the plot above, computeMatrix can be run with either one of the two modes: __scaled regions__ or __reference point__.

Please see the [example figures](#examples) down below for explanations of parameters and options.


### Output files
  + __obligatory__: zipped matrix of values to be used with _heatmapper_ and/or _profiler_
  + __optional__  (can also be generated with heatmapper or profiler in case you forgot to produce them in the beginning):
    - BED-file of the regions sorted according to the calculated values
    - list of average values per genomic bin
    - matrix of values per genomic bin per genomic interval


<a name="heatmapper"/>
## heatmapper
The heatmapper depicts values extracted from the [bigWig][] file for each genomic region individually.
It requires the output from computeMatrix and most of its options are related to tweeking the visualization only. The values calculated by computeMatrix are not changed.

Definitely check the example at the bottom of the page to get a feeling for how many things you can tune.


<a name="profiler"/>
## profiler
This tool plots the average enrichments over all genomic regions supplied to computeMarix. It is a very useful complement to the heatmapper, especially in cases when you want to compare the scores for many different groups. Like heatmapper, profiler does not change the values that were compute by computeMatrix, but you can choose between many different ways to color and display the plots.


<a name="examples"/>
## Example figures

Here you see a typical, not too pretty example of a heatmap. We will use this example to explain several features of computeMatrix and heatmapper, so do take a closer look.

### Heatmap with all genes scaled to the one size and user-specified groups of genes
![Heatmap](https://raw.github.com/fidelram/deepTools/master/examples/visual_hm_DmelPolII.png "Heatmap of RNA Polymerase II ChIP-seq")

The plot was produced with the following commands:

    $ /deepTools-1.5/bin/computeMatrix scale-regions -R Dm.genes.indChromLabeled.bed -S PolII.bw --beforeRegionStartLength 500 --afterRegionStartLength 500 --regionBodyLength 1500 --binSize 10 -p20 --outFileName PolII_matrix_scaledGenes --sortRegions no
    $ /deepTools-1.5/bin/heatmapper -m PolII_matrix_scaledGenes --outFileName PolII_indChr_scaledGenes.pdf -T "Pol II" --whatToShow "heatmap only"

As you can see, all genes have been scaled to the same size and the (mean) values per bin size (10 bp) are colored accordingly. In addition to the gene bodies, we added 500 bp up- and down-stream of the genes.

This is what you would have to select to achieve the same result within Galaxy:

###### computeMatrix
![computeMatrixGal01](https://raw.github.com/fidelram/deepTools/master/examples/visual_computeMatrix03.png "deepTools Galaxy screenshot of computeMatrix")
![computeMatrixGal02](https://raw.github.com/fidelram/deepTools/master/examples/visual_computeMatrix02.png "deepTools Galaxy screenshot of computeMatrix (advanced options)")

###### heatmapper
![computeMatrixGal03](https://raw.github.com/fidelram/deepTools/master/examples/visual_heatmapper.png "deepTools Galaxy screenshot of heatmapper (advanced options)")

###### main difference between computeMatrix usage on the command line and Galaxy: the input of the regions file (BED) 

Note that we supplied just _one_ BED-file via the command line whereas in Galaxy we indicated three different files (one per chromosome).

On the command line, the program expects a BED file where different groups of genomic regions are concatenated into one file, where the beginning of each group should be indicated by "#group name".
The BED-file that was used here, contained 3 such lines and could be prepared as follows:
    
     $ grep ^chr2 AllGenes.bed > Dm.genes.indChromLabeled.bed

     $ echo "#chr2" >> Dm.genes.indChromLabeled.bed
     
     $ grep ^chr3 AllGenes.bed >> Dm.genes.indChromLabeled.bed
     
     $ echo "#chr3" >> Dm.genes.indChromLabeled.bed
     
     $ grep ^chrX AllGenes.bed >> Dm.genes.indChromLabeled.bed
     
     $ echo "#chrX" >> Dm.genes.indChromLabeled.bed
     
In Galaxy, you can simply generate three different data sets starting from a whole genome list by using the "Filter" tool three times:
1. c1=="chr2" --> Dm.genes.chr2.bed
2. c1=="chr3" --> Dm.genes.chr3.bed
3. c1=="chrX" --> Dm.genes.chrX.bed
4. 

#### Important parameters for optimizing the visualization
1. __sorting of the regions__: The default of heatmapper is to sort the values descendingly. You can change that to ascending, no sorting at all or according to the size of the region (Using the `--sort` option on the command line or advanced options in Galaxy). We strongly recommend to leave the sorting option at "no sorting" for the intitial computeMatrix step.
2. __coloring__: The default coloring by heatmapper is done using the python colormap "RdYlBu", but this can be changed (--colorMap on the command line, advanced options within Galaxy).
4. __dealing with missing data__: You have certainly noticed that some gene bodies are depicted as white lines within the otherwise colorful mass of genes. Those regions are due to genes that, for whatever reason, did not have any read coverage in the bigWig file. There are several ways to handle these cases:
    + __--skipZeros__ this is useful when your data actually has a quite nice coverage, but there are 2 or 3 regions where you deliberately filtered out reads or you don't expect any coverage (e.g. hardly mappable regions). This will only work if the entire region does not contain a single value. 
    + __--missingDataAsZero__ this option allows computeMatrix do interpret missing data points as zeroes. Be aware of the changes to the average values that this might cause.
    + __--missingDataColor__ this is in case you have very sparce data or were missing values make sense (e.g. when plotting methylated CpGs - half the genome should have no value). This option then allows you to pick out your favorite color for those regions. The default is black (was white when the above shown image was produced).


### Summary plots
Here's the __profiler__ plot corresponding to the heatmap above. There's one major difference though - do you spot it?

![Profile](https://raw.github.com/fidelram/deepTools/master/examples/visual_profiler_DmelPolII.png "Meta-gene profile of RNA Polymerase II")

We used the same [BED][] file(s) as for the heatmap, hence the 3 different groups (1 per chromosome). However, this time we used computeMatrix not with _scale-regions_ but with _reference-point_ mode.

    $ /deepTools-1.5/bin/computeMatrix reference-point --referencePoint TSS -R Dm.genes.indChromLabeled.bed -S PolII.bw -b 1000 -a 1000 -bs 10 -p20 --outFileName PolII_matrix_indChr_refPoint --missingDataAsZero --sortRegions no
    $ /deepTools-1.5/bin/profiler --matrixFile PolII_matrix_indChr_refPoint --outFileName profile_PolII_indChr_refPoint.pdf --plotType fill --startLabel "TSS" -T "Pol II around TSS" --yAxisLabel "mean Pol II coverage" --onePlotPerGroup
 
When you compare the profiler commands with the heatmapper commands, you also notice that we made use of many more labelling options here, e.g. `--yAxisLabel` and a more specific title via `-T`


This is how you would have obtained this plot in Galaxy (only the part that's _different_ from the above shown command for the scale-regions version is shown):

###### computeMatrix
![computeMatrixGal04](https://raw.github.com/fidelram/deepTools/master/examples/visual_computeMatrix03.png "deepTools Galaxy screenshot of computeMatrix for profiles in reference-point mode")
###### profiler
![computeMatrixGal04](https://raw.github.com/fidelram/deepTools/master/examples/visual_profiler_Gal.png "deepTools Galaxy screenshot of profilers in reference-point mode")



### Heatmap with all genes scaled to the one size and kmeans clustering

Instead of supplying groups of regions on your own, you can use the clustering function of heatmapper to get a first impression whether the signal of your experiment can be easily clustered into two or more groups of similar signal distribution.

Have  a look at this example with two clusters:

![kmeans](https://raw.github.com/fidelram/deepTools/master/examples/heatmaps_kmeans_Pol_II.png "Heatmap of RNA Polymerase II ChIP-seq divided into two clusters.")

The plot was produced with the following commands:

    $ /deepTools-1.5/bin/computeMatrix reference-point -R Dm.genes.indChromLabeled.bed -S PolII.bw --beforeRegionStartLength 500 --afterRegionStartLength 500o  --binSize 50 -p20 --outFileName PolII_matrix_TSS
    $ /deepTools-1.5/bin/heatmapper -m PolII_matrix_TSS --kmeans 2 --outFileName PolII_two_clusters.pdf -T "Pol II" --sortUsing region_length --whatToShow "heatmap only"

When the `--kmeans` option is chosen and more than 0 clusters are specified, heatmapper will run the [k-means][] clustering algorithm. In this example _Drosophila m._ genes were divided into two clusters separating those genes with Pol II at the promoter region (top) from those genes without Poll II at the promoter (bottom).
Please note that the clustering will only work if the initial BED-file used with computeMatrix contained only _one_ group of genes (i.e. all genes, without any hash tags separating them)

The genes belonging to each cluster can be obtained by via `--outFileSortedRegions` on the command line and "advanced output options in Galaxy". On the command line, this will result in a BED file where the groups are separated by a hash tag. In Galaxy, you will obtain individual data sets per cluster.

To have a better control on the clustering it is recommended to load the matrix raw data into __specialized software like [cluster3] or [R]__. You can obtain the matrix via the option `--outFileNameMatrix` on the command line and by the "advanced output options" in Galaxy. The order of the rows is the same as in the ouput of the `--outFileSortedRegions` BED file.


-----------------------------------------------------------------------------------
[BAM]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "binary version of a SAM file; contains all information about aligned reads"
[BED]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "simple text file of genomic regions (chr, start, end)"
[SAM]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file containing all information about aligned reads"
[bigWig]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "binary version of a bedGraph file; contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[bedGraph]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file that contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[FASTQ]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file of raw reads (almost straight out of the sequencer)"
### References
[Benjamini and Speed]: http://nar.oxfordjournals.org/content/40/10/e72 "Nucleic Acids Research (2012)"
[Diaz et al.]: http://www.degruyter.com/view/j/sagmb.2012.11.issue-3/1544-6115.1750/1544-6115.1750.xml "Stat. Appl. Gen. Mol. Biol. (2012)"
[k-means]: http://en.wikipedia.org/wiki/K-means_clustering
[cluster3]: http://bonsai.hgc.jp/~mdehoon/software/cluster/
[R]: http://www.r-project.org/

This tool is developed by the [Bioinformatics Facility](http://www1.ie-freiburg.mpg.de/bioinformaticsfac) at the [Max Planck Institute for Immunobiology and Epigenetics, Freiburg](http://www1.ie-freiburg.mpg.de/).
