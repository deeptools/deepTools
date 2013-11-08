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
Here you see a typical, not too pretty example of a heatmap. We will use this example to explain several features of heatmapper, so do take a closer look.

![Heatmap](https://raw.github.com/fidelram/deepTools/master/examples/visual_hm_DmelPolII.png "Heatmap of RNA Polymerase II ChIP-seq")

The plot was produced with the following commands:

    $ /deepTools-1.5/bin/computeMatrix scale-regions -R Dm.genes.indChromLabeled.bed -S PolII.bw --beforeRegionStartLength 500 --afterRegionStartLength 500 --regionBodyLength 1500 --binSize 10 -p20 --outFileName PolII_matrix_scaledGenes --sortRegions no
    $ /deepTools-1.5/bin/heatmapper -m PolII_matrix_scaledGenes --outFileName PolII_indChr_scaledGenes.pdf -T "Pol II" --whatToShow "heatmap only"


As you can see, all genes have been scaled to the same size and the (mean) values per bin size (10 bp) are colored accordingly. In addition to the gene bodies, we added 500 bp up- and down-stream of the genes.

1. __sorting of the regions__ The default of heatmapper is to sort the values descendingly. Using the --sort option, you can change that to ascending, no sorting at all or according to the size of the region.
2. __coloring__ The default coloring by heatmapper is done using the python colormap "RdYlBu".
4. __stretches of white__ You certainly have noticed that some gene bodies are depicted as white lines within the otherwise colorful mass of genes. Those regions are due to genes that, for whatever reason, did not have any read coverage scores in the bigWig file. There are several ways to handle these cases:
    + __--skipZeros__ this is useful when your data actually has a quite nice coverage, but there are 2 or 3 regions where you deliberately filtered out reads or you don't expect any coverage (e.g. hardly mappable regions). This will only work if the entire region does not contain a single value. 
    + __--missingDataAsZero__ this option allows computeMatrix do interpret missing data points as zeroes. Be aware of the changes to the average values that this might cause.
    + __--missingDataColor__ this is in case you have very sparce data or were missing values make sense (e.g. when plotting methylated CpGs - half the genome should have no value). This option then allows you to pick out your favorite color for those regions. The default is white.
3. __1 group per chromosome__ To tell computeMatrix that there are different groups of genomic regions, the [BED][]-file needs to contain the name of the group, preceded by a hash tag __at the end of each group__. The BED-file that was used here, contained 3 such lines and could be prepared as follows:
    
     $ grep ^chr2 AllGenes.bed > Dm.genes.indChromLabeled.bed

     $ echo "#chr2" >> Dm.genes.indChromLabeled.bed
     
     $ grep ^chr3 AllGenes.bed >> Dm.genes.indChromLabeled.bed
     
     $ echo "#chr3" >> Dm.genes.indChromLabeled.bed
     
     $ grep ^chrX AllGenes.bed >> Dm.genes.indChromLabeled.bed
     
     $ echo "#chrX" >> Dm.genes.indChromLabeled.bed


Here's the profiler plot corresponding to the heatmap above. There's one major difference though - do you spot it?

![Profile](https://raw.github.com/fidelram/deepTools/master/examples/visual_profiler_DmelPolII.png "Meta-gene profile of Rna Polymerase II")

We used the same [BED][] file as for the heatmap, hence the 3 different groups (1 per chromosome). The major difference was that we used computeMatrix not with _scale-regions_ but with _reference-point_ mode.

    $ /deepTools-1.5/bin/computeMatrix reference-point --referencePoint TSS -R Dm.genes.indChromLabeled.bed -S PolII.bw -b 1000 -a 1000 -bs 10 -p20 --outFileName PolII_matrix_indChr_refPoint --missingDataAsZero --sortRegions no
    $ /deepTools-1.5/bin/profiler --matrixFile PolII_matrix_indChr_refPoint --outFileName profile_PolII_indChr_refPoint.pdf --plotType fill --startLabel "TSS" -T "Pol II around TSS" --yAxisLabel "mean Pol II coverage" --onePlotPerGroup
 
When you compare the profiler commands with the heatmapper commands, you also notice that we made use of many more labelling options here, e.g. --yAxisLabel and a more specific title via -T



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


This tool is developed by the [Bioinformatics Facility](http://www1.ie-freiburg.mpg.de/bioinformaticsfac) at the [Max Planck Institute for Immunobiology and Epigenetics, Freiburg](http://www1.ie-freiburg.mpg.de/).
