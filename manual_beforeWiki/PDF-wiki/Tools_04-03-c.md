
When the `--kmeans` option is chosen and more than 0 clusters are specified, heatmapper will run the [k-means][] clustering algorithm. In this example, we wanted to divide _Drosophila melanogaster_ genes into two clusters. As you can see above, the algorithm nicely identified two groups - one with mostly those genes with lots of Pol II at the promoter region (top) from those genes without Poll II at the promoter (bottom).
Please note that the clustering will only work if the initial BED-file used with computeMatrix contained only _one_ group of genes.

The genes belonging to each cluster can be obtained by via `--outFileSortedRegions` on the command line and "advanced output options in Galaxy". On the command line, this will result in a BED file where the groups are separated by a hash tag. In Galaxy, you will obtain individual data sets per cluster.

To have a better control on the clustering it is recommended to load the matrix raw data into __specialized software like [cluster3] or [R]__. You can obtain the matrix via the option `--outFileNameMatrix` on the command line and by the "advanced output options" in Galaxy. The order of the rows is the same as in the output of the `--outFileSortedRegions` BED file.

[Benjamini and Speed]: http://nar.oxfordjournals.org/content/40/10/e72 "Nucleic Acids Research (2012)"
[Diaz et al.]: http://www.degruyter.com/view/j/sagmb.2012.11.issue-3/1544-6115.1750/1544-6115.1750.xml "Stat. Appl. Gen. Mol. Biol. (2012)"
[k-means]: http://en.wikipedia.org/wiki/K-means_clustering
[cluster3]: http://bonsai.hgc.jp/~mdehoon/software/cluster/
[R]: http://www.r-project.org/
[IGV]: http://www.broadinstitute.org/igv/ "Integrative Genome Browser developed by the Broad Institute"
[bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml