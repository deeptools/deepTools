Frequently asked questions
===========================

#### [How can I do...?](#HowTo)
* [I have downloaded/received a BAM file - how do I generate a file I can look at in a Genome Browser?](#FASTQ2IGV)
* [How can I assess the reproducibility of my sequencing replicates?](#repCorr)
* [How do I get an input-normalized ChIP-seq coverage file?](#InputNorm)
* [How can I compare the ChIP strength for different ChIP experiments?](#fprint)
* [How do I know whether my sample is GC biased? And if yes, how do I correct for it?](#GC)
* [How do I get a (clustered) heatmap of sequencing-depth-normalized read coverages around the transcription start site of all genes?](#HM)
* [How can I compare the average signal for X- and autosomal genes for 2 or more different sequencing experiments](#profiler)

#### [Galaxy-specific questions](#GalSpecific)
* [I've reached my quota - what can I do to save some space?](#quota)
* [What is the best way to integrate the deepTools results with other downstream analyses (outside of Galaxy)?](#integrate)
* [How can I determine basic parameters of a BAM file?](#BAMparams)
* [How can I get a measure of mappability and GC content for my genome of interest?](#MapGC)

#### [General deepTools-related questions](#general)
* [I just want to try out a tool, how can I optimize the computation time?](#compTime)
* [When should I exclude regions from computeGCbias?](#excludeGC)
* [Does it speed up the computation if I limit bamCorrelate to one chromosome, but keep the same numbers and sizes of sampling bins?](#bamCorrelateLimit)
* [How do I calculate the effective genome size for an organism that's not in your list?](#effGenomeSize)
* [How can I increase the resolution of the heatmap?](#hmresolution)


###### [Back to deepTools Galaxy](https://deeptools.ie-freiburg.mpg.de)

###### [Go to the general help page](https://github.com/fidelram/deepTools/blob/master/manual/GalaxyHelp.md)

----------------------------------------------------------------------------------------------------


How can I do...? <a name="HowTo"></a>
--------------------------------------------------
This section is meant to give you quick guidance to specific tasks you may want to perform. There are many more ways in which you can use [deepTools Galaxy][] than those described here, so be creative once you're comfortable with using them. For detailed explanations of what the tools do, follow the links.
All of these tipps assume that you're working on the [deepTools Galaxy][] and have uploaded your files.

#### I have downloaded/received a [BAM][] file - how do I generate a file I can look at in a Genome Browser?<a name="FASTQ2IGV"></a>
* tool: [bamCoverage][]
* input: your BAM file

![GalHow_bamCoverage](https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamCoverage.png "deepTools Galaxy screenshot of bamCoverage usage and output")
--------------------------------------------------

#### How can I assess the reproducibility of my sequencing replicates?<a name="repCorr"></a>
* tool: [bamCorrelate][]
* input: BAM files
    * you can compare as many samples as you want - the more you put at the same time, the longer the computation takes
* output: heatmap of correlations - the closer two samples are to each other, the more similar their read coverages

![GalHow_bamCorrelate](https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamCorrelate.png "deepTools Galaxy screenshot of bamCorrelate usage and output")
--------------------------------------------------

#### How do I get an input-normalized ChIP-seq coverage file?<a name="InputNorm"></a>

1. you need two BAM files: one for the input, one for the ChIP-seq experiment
2. use the tool [bamCompare][] with ChIP = treatment, input = control sample

![GalHow_bamCompare](https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamCompare.png "deepTools Galaxy screenshot of bamCompare usage and output")
--------------------------------------------------

#### How can I compare the ChIP strength for different ChIP experiments?<a name="fprint"></a>
* tool: [bamFingerprint][]
* input: as many BAM files as you'd like to compare. Make sure you get all the labels right!

![GalHow_Fingerprint](https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamFingerprint.png "deepTools Galaxy screenshot of bamFingerprint")
--------------------------------------------------

#### How do I know whether my sample is GC biased? And if yes, how do I correct for it?<a name="GC"></a>
* you need a BAM file of your sample in question
* use the tool [computeGCbias][] on that BAM file (default settings, just make sure your reference genome and genome size are matching)

![GalHow_computeGC](https://raw.github.com/fidelram/deepTools/master/examples/GalHow_computeGCbias.png "deepTools Galaxy screenshot of computeGCbias")

* have a look at the image that is produced and compare it to the examples [here](https://github.com/fidelram/deepTools/blob/master/manual/QC.md#computeGCbias)
* if your sample shows an almost linear increase in exp/obs coverage (on the log scale of the lower plot), then you should consider correcting the GC bias - _if_ you think that the biological interpretation of this data would otherwise be compromised (e.g. by comparing it to another sample that does not have an inherent GC bias)
    
    + the GC bias can be corrected with the tool [correctGCbias][] using the second output of the computeGCbias tool that you had to run anyway
    + CAUTION!! correctGCbias will add reads to otherwise depleted regions (typically GC-poor regions), that means that you should __not__ remove duplicates in any downstream analyses based on the GC-corrected BAM file (we therefore recommend to remove duplicates before doing the correction so that only those duplicate reads are kept that were produced by the GC correction procedure)

![GalHow_correctGC](https://raw.github.com/fidelram/deepTools/master/examples/GalHow_correctGCbias.png "deepTools Galaxy screenshot of correctGCbias usage and output")
--------------------------------------------------

#### How do I get a (clustered) heatmap of sequencing-depth-normalized read coverages around the transcription start site of all genes?<a name="HM"></a>
* if you want to start with a BAM file, begin by _generating the normalized read coverages_ using the tool [bamCoverage][] with the option "normalize to 1x sequencing depth" (make sure that you indicate the correct genome size) (1)
* you also need a BED or INTERVAL file of genes (you can obtain one via "Get Data" --> "UCSC main table browser" --> group: "Genes and Gene Predictions" --> (e.g.) "RefSeqGenes" --> send to Galaxy (2)


![GalHow_clustHM01](https://raw.github.com/fidelram/deepTools/master/examples/GalHow_clustHM01.png "deepTools Galaxy screenshot of how to get a list of genes from UCSC")

* use [computeMatrix][] with the coverage file generated in (1) and the BED file from (2), indicate "reference-point" and whatever other option you would like to tune (3)

![GalHow_clustHM02](https://raw.github.com/fidelram/deepTools/master/examples/GalHow_clustHM02.png "deepTools Galaxy screenshot of computeMatrix for profiles in reference-point mode with output")

* use the output from (3) with [heatmapper][] (if you would like to cluster the signals, choose "kmeans clustering" (last option of "advanced options") with a reasonable number of clusters)

![GalHow_clustHM03](https://raw.github.com/fidelram/deepTools/master/examples/GalHow_clustHM03.png "deepTools Galaxy screenshot of heatmapper usage and output")
--------------------------------------------------

#### How can I compare the average signal for X- and autosomal genes for 2 or more different sequencing experiments?<a name="profiler"></a>
* you need two __BED files__: one with X-chromosomal and one with autosomal genes (1)
    * you can download a full list of genes via "Get Data" --> "UCSC main table browser" --> group:"Genes and Gene Predictions" --> tracks: (e.g.) "RefSeqGenes" --> send to Galaxy
    * then filter the full list twice using the tool "Filter data on any column using simple expressions" 
        - first use the expression: c1=="chrX" to filter the list of all genes --> this will generate a list of X-linked genes
        - then re-run the filtering, now with c1!="chrX" which will generate a list of genes that do not belong to chromosome X (!= indicates "not matching")
* you need __bigWig files__ for each experiment (in case you only have BAM files, run [bamCoverage][] on every BAM file first) (2)
* use [computeMatrix][] for each signal file (bigWig) (you only need to specify all the parameters once, then use the re-run button underneath the first data set and just replace the signal file with the next one) (3)
    * supply both filtered BED files (click on "Add new regions to plot" once) and label them
    * indicate the corresponding signal file
    * make sure to __re-name__ every data set in the history once computeMatrix is done so that you can easily keep track of which matrix was based on which bigWig file (you can always find these information by clicking on the i-button in the respective data set)
* now use [profiler][] for every file you generated with computeMatrix (4)
    * important: display the "advanced output options" and select "save the data underlying the average profile" --> this will generate a table in addition to the summary plot images
* now you have at least 2 separate images of profiles - one for each bigWig file - you can either leave it like this or use another script that will plot all the summary plots in one image at once (5)
    * this tool is called "Plotting the summary plots for multiple signals"
    * it uses the tables generated by profiler in (4)
    * for each group of genes (in this case, X and autosomal genes = 2 groups), you can assign a color
--------------------------------------------------



Galaxy-specific questions <a name="GalSpecific">
--------------------------------------------------

#### I've reached my quota - what can I do to save some space? <a name="quota"></a>
1. make sure that all the data sets you deleted are __permanently__ eliminated from our disks: go to the history option button and select "Purge deleted data sets", then hit the "refresh" button on top of your history panel
2. download all data sets for which you've completed the analysis, then remove the data sets (click on the "x" and then make sure they're purged (see above)


#### What's the best way to integrate the deepTools results with other downstream analyses (outside of Galaxy) <a name="integrate"></a>
* you can save all the data tables underlying every image produced by deepTools, i.e. if you would like to plot the average profiles in a different way, you could download the corresponding data (after ticking the profiler option at "advanced output options") and import them into R, Excel, GraphPadPrism etc.


#### How can I determine basic parameters of a BAM file, such as the number of reads, read length, duplication rate and average DNA fragment length? <a name="BAMparams"></a>

Eventhough [MACS][] is meant to do peak calling for you, it also outputs a number of useful information such as those listed above.
Simply run MACS on the BAM file that you would like to gain the information for and check the .xls file from the MACS output. It will list:

* tag length = read length
* duplication rate
* number of tags = number of reads
* d = distance = average DNA fragment size


#### How can I get a measure of mappability and GC content for my genome of interest?<a name="MapGC"></a>
In the tools panel, go to "Get Data" --> "UCSC main table browser" --> select your organism and the genome assembly you're interested in --> select "Mapping and Sequencing Tracks" (group menu) --> select the track



General deepTools-related questions <a name="general">
--------------------------------------------------------------

#### How can I test a tool with little computation time? <a name="compTime"></a>
* when you're playing around with the tools to see what kinds of results they will produce: limit the operation to one chromosome only to __save computation time__! ("advanced output options" --> "Region of the genome to limit the operation to")


#### When should I exclude regions from computeGCbias? <a name="excludeGC"></a>
In general, we recommend that you should only correct for GC bias (using computeGCbias followed by correctGCbias) if you observe that the majority of the genome (the region between 30-60%) is continuously GC-biased __and__ you want to compare this sample with another sample that is not GC-biased.

Sometimes, a certain GC bias is expected, for example for ChIP samples of H3K4me3 in mammalian samples where GC-rich promoters are expected to be enriched. To not confound the GC bias caused by the library preparation with the inherent, expected GC bias, we incorporated the possibility to supply a file of regions to computeGCbias that will be excluded from the GC bias calculation. This file should typically contain those regions that one expects to be significantly enriched per se. This way, the computeGCbias will focus on background regions.


#### Does it speed up the computation if I limit bamCorrelate to one chromosome, but keep the same numbers and sizes of sampling bins?<a name="bamCorrelateLimit"></a>
Yes. However, the way bamCorrelate (and all the other deepTools handle the option "limit the computation to a specific region" is as follows: first, the _entire_ genome represented in the BAM file will be regarded and sampled, _then_ all the regions or sampled bins that do not overlap with the region indicated by the user will be discarded. This means that if you wanted 10,000 bins to be sampled and you focus on, let's say, chromosome 2, the final computation will not be performed on the whole set of 10,000 bins, but only on those bins that overlap with chromosome 2.


#### How do I calculate the effective genome size for an organism that's not in your list?<a name="effGenomeSize"></a>
This is something you will have to find a solution outside of deepTools at the moment. We suggest to run faCount from UCSC tools. If you used multi-read alignment (e.g. with bowtie2), then you can use that tool to report the total number of bases as well as the number of unmapped bp, indicated by 'N'. The effective genome size is the total number of reads minus the number of 'N'.


#### The heatmap I generated looks very "coarse", I would like a much more fine-grained image. <a name="hmresolution"></a>
* decrease the __bin size__ when generating the matrix using computeMatrix
  * go to "advanced options" --> "Length, in base pairs, of the non-overlapping bin for averaging the score over the regions length" --> define a smaller value, e.g. 50 or 25 bp
* make sure, however, that you used a sufficiently small bin size when calculating the bigWig file, though (if generated with deepTools, you can check the option "bin size")

----------------------------------------------------------------
##### [Back to general deepTools Galaxy help page](https://github.com/fidelram/deepTools/blob/master/manual/GalaxyHelp.md#deepTools)
----------------------------------------------------------------
[Galaxy]: http://galaxyproject.org/ "General Galaxy platform from Penn State"
[GEO]: http://www.ncbi.nlm.nih.gov/geo/ "GEO database"
[Roadmap project]: http://www.roadmapepigenomics.org/data "Roadmap web site"
[UCSC]: http://genome.ucsc.edu/ "UCSC Genome web site"
[BioMart]: http://www.biomart.org/ "Biomart web site"
[deepTools Galaxy]: http://deeptools.ie-freiburg.mpg.de/ "deepTools Galaxy at the Max-Planck-Institute of Immunobiology and Epigenetics"

[bamCorrelate]: https://github.com/fidelram/deepTools/blob/master/manual/QC.md
[bamFingerprint]: https://github.com/fidelram/deepTools/blob/master/manual/QC.md
[computeGCBias]: https://github.com/fidelram/deepTools/blob/master/manual/QC.md
[bamCoverage]: https://github.com/fidelram/deepTools/blob/master/manual/normalizations.md
[bamCompare]: https://github.com/fidelram/deepTools/blob/master/manual/normalizations.md
[correctGCbias]: https://github.com/fidelram/deepTools/blob/master/manual/normalizations.md
[computeMatrix]: https://github.com/fidelram/deepTools/blob/master/manual/visualizations.md
[heatmapper]: https://github.com/fidelram/deepTools/blob/master/manual/visualizations.md
[profiler]: https://github.com/fidelram/deepTools/blob/master/manual/visualizations.md
[MACS]: http://www.ncbi.nlm.nih.gov/pubmed/22936215 "How to use MACS, Nature Protocols"
[CCAT]: http://www.ncbi.nlm.nih.gov/pubmed/20371496 "CCAT original publication"
[SICER]: http://bioinformatics.oxfordjournals.org/content/25/15/1952.full "SICER original publication"

[BAM]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "binary version of a SAM file; contains all information about aligned reads"
[SAM]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file containing all information about aligned reads"
[bigWig]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "binary version of a bedGraph file; contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[bedGraph]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file that contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[FASTQ]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file of raw reads (almost straight out of the sequencer)"

 
