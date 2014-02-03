### 2nd example: Summary plots with all genes scaled to the one size and user-specified groups of genes

Here's the __profiler__ plot corresponding to the heatmap above. There's one major difference though - do you spot it?

![Profile](https://raw.github.com/fidelram/deepTools/master/examples/visual_profiler_DmelPolII.png "Meta-gene profile of RNA Polymerase II")

We used the same BED file(s) as for the heatmap, hence the 3 different groups (1 per chromosome). However, this time we used computeMatrix not with _scale-regions_ but with _reference-point_ mode.

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
 
When you compare the profiler commands with the heatmapper commands, you also notice that we made use of many more labeling options here, e.g. `--yAxisLabel` and a more specific title via `-T`


This is how you would have obtained this plot in Galaxy (only the part that's _different_ from the above shown command for the scale-regions version is shown):

##### computeMatrix

<a href="https://raw.github.com/fidelram/deepTools/master/examples/visual_computeMatrix03.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/visual_computeMatrix03.png" Title="deepTools Galaxy screenshot of computeMatrix for profiles in reference-point mode" />
</a>

##### profiler
<a href="https://raw.github.com/fidelram/deepTools/master/examples/visual_profiler_Gal.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/visual_profiler_Gal.png" Title="deepTools Galaxy screenshot of profiler in reference-point mode" />
</a>


[Benjamini and Speed]: http://nar.oxfordjournals.org/content/40/10/e72 "Nucleic Acids Research (2012)"
[Diaz et al.]: http://www.degruyter.com/view/j/sagmb.2012.11.issue-3/1544-6115.1750/1544-6115.1750.xml "Stat. Appl. Gen. Mol. Biol. (2012)"
[k-means]: http://en.wikipedia.org/wiki/K-means_clustering
[cluster3]: http://bonsai.hgc.jp/~mdehoon/software/cluster/
[R]: http://www.r-project.org/
[IGV]: http://www.broadinstitute.org/igv/ "Integrative Genome Browser developed by the Broad Institute"
[bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml