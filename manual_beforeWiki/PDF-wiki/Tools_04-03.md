### 3rd example: Heatmap with all genes scaled to the one size and kmeans clustering

Instead of supplying groups of regions on your own, you can use the clustering function of heatmapper to get a first impression whether the signal of your experiment can be easily clustered into two or more groups of similar signal distribution.

Have  a look at this example with two clusters. The values correspond to log2ratios(ChIP/input) from a ChIP-seq experiment for RNA Polymerase II in _Drosophila melanogaster_:

<a href="https://raw.github.com/fidelram/deepTools/master/examples/heatmaps_kmeans_Pol_II_small.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/heatmaps_kmeans_Pol_II_small.png" Title="Heatmap of RNA Polymerase II ChIP-seq divided into two clusters." />
</a>

The plot was produced with the following commands:

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

In Galaxy, these are the screenshots from the commands for computeMatrix and heatmapper:

[Benjamini and Speed]: http://nar.oxfordjournals.org/content/40/10/e72 "Nucleic Acids Research (2012)"
[Diaz et al.]: http://www.degruyter.com/view/j/sagmb.2012.11.issue-3/1544-6115.1750/1544-6115.1750.xml "Stat. Appl. Gen. Mol. Biol. (2012)"
[k-means]: http://en.wikipedia.org/wiki/K-means_clustering
[cluster3]: http://bonsai.hgc.jp/~mdehoon/software/cluster/
[R]: http://www.r-project.org/
[IGV]: http://www.broadinstitute.org/igv/ "Integrative Genome Browser developed by the Broad Institute"
[bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml