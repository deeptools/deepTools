Step-by-step protocols
========================

.. contents:: 
    :local:


How can I do...? <a name="HowTo"></a>
--------------------------------------------------

This section should give you a quick overview of how to do many common tasks. We're using screenshots from Galaxy here, so if you're using the command-line version then you can easily follow the given examples by typing the program name and the help option (e.g. /deepTools/bin/bamCoverage --help), which will show you all the parameters and options (most of them named very similarly to those in Galaxy).

For each "recipe" here, you will find the screenshot of the tool and the input parameters on the left hand side (we marked non-default, _user-specified entries_) and screenshots of the output on the right hand side. Do let us know if you spot things that are missing, should be explained better, or are simply confusing!

There are many more ways in which you can use [deepTools Galaxy][] than those described here, so be creative once you're comfortable with using them. For detailed explanations of what the tools do, follow the links.

> All recipes assume that you have uploaded your files into a Galaxy instance with a deepTools installation, e.g. [deepTools Galaxy][].

> If you would like to try out the protocols with __sample data__, go to [deepTools Galaxy][]  &rarr; "Shared Data"  &rarr; "Data Libraries"  &rarr; "deepTools Test Files". Simply select BED/BAM/bigWig files and click, "to History". You can also download the test datasets by clicking "Download" at the top.

___________________________________
<a name="FASTQ2IGV"></a>
#### I have downloaded/received a [BAM][] file - how do I generate a file I can look at in a genome browser?

* tool: [bamCoverage][]
* input: your BAM file

Note: BAM files can also be viewed in genome browsers, however, they're large and tend to freeze the applications. Generating bigWig files of read coverages will help you a lot in this regard. In addition, if you have more than one sample you'd like to look at, it is helpful to normalize all of them to 1x sequencing depth.

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamCoverage.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamCoverage.png" Title="deepTools Galaxy screenshot of bamCoverage usage and output" />
</a>

___________________________________
<a name="repCorr"></a>
#### How can I assess the reproducibility of my sequencing replicates?

* tool: [multiBamCoverage][]
* input: BAM files
    * you can compare as many samples as you want, though the more you use the longer the computation will take

* output: heatmap of correlations - the closer two samples are to each other, the more similar their read coverages will be

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_multiBamCoverage.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_multiBamCoverage.png" Title="deepTools Galaxy screenshot of multiBamCoverage usage and output" />
</a>

___________________________________
<a name="GC"></a>
#### How do I know whether my sample is GC biased? And if it is, how do I correct for it?

* you need a BAM file of your sample
* use the tool [computeGCbias][] on that BAM file (default settings, just make sure your reference genome and genome size are matching)

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_computeGCbias.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_computeGCbias.png" Title="deepTools Galaxy screenshot of computeGCbias" />
</a>


* have a look at the image that is produced and compare it to the examples [here](https://github.com/fidelram/deepTools/wiki/QC)
* if your sample shows an almost linear increase in exp/obs coverage (on the log scale of the lower plot), then you should consider correcting the GC bias - _if_ you think that the biological interpretation of this data would otherwise be compromised (e.g. by comparing it to another sample that does not have an inherent GC bias)

    + the GC bias can be corrected with the tool [correctGCbias][] using the second output of the computeGCbias tool that you had to run anyway

    + CAUTION!! correctGCbias will add reads to otherwise depleted regions (typically GC-poor regions), that means that you should __not__ remove duplicates in any downstream analyses based on the GC-corrected BAM file (we therefore recommend removing duplicates before doing the correction so that only those duplicate reads are kept that were produced by the GC correction procedure)

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_correctGCbias.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_correctGCbias.png" Title="deepTools Galaxy screenshot of correctGCbias usage and output" />
</a>

___________________________________
<a name="InputNorm"></a>
#### How do I get an input-normalized ChIP-seq coverage file?

* input: you need two BAM files, one for the input and one for the ChIP-seq experiment
* tool: [bamCompare][] with ChIP = treatment, input = control sample

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamCompare.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamCompare.png" Title="deepTools Galaxy screenshot of bamCompare usage and output" />
</a>

___________________________________
<a name="fprint"></a>
#### How can I compare the ChIP strength for different ChIP experiments?

* tool: [plotFingerprint][]
* input: as many BAM files as you'd like to compare. Make sure you get all the labels right!

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_plotFingerprint.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_plotFingerprint.png" Title="deepTools Galaxy screenshot of plotFingerprint" />
</a>

___________________________________
<a name="HM"></a>
#### How do I get a (clustered) heatmap of sequencing-depth-normalized read coverages around the transcription start site of all genes?

* tools: [computeMatrix][], then [heatmapper][]
* inputs:
    * 1 bigWig file of normalized read coverages (e.g. the result of bamCoverage or bamCompare)
    * 1 BED or INTERVAL file of genes, e.g. obtained through Galaxy via "Get Data" &rarr; "UCSC main table browser" &rarr; group: "Genes and Gene Predictions" &rarr; (e.g.) "RefSeqGenes" &rarr; send to Galaxy (see screenshots below)

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_clustHM01.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_clustHM01.png" Title="deepTools Galaxy screenshot of how to get a list of genes from UCSC" />
</a>

* use [computeMatrix][] with the bigWig file and the BED file
* indicate "reference-point"  (and whatever other option you would like to tune, see screenshot below)

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_clustHM02.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_clustHM02.png" Title="deepTools Galaxy screenshot of computeMatrix for profiles in reference-point mode with output" />
</a>


* use the output from computeMatrix with [heatmapper][]
    * if you would like to cluster the signals, choose "k-means clustering" (last option of "advanced options") with a reasonable number of clusters (usually between 2 to 7)

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_clustHM03.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_clustHM03.png" Title="deepTools Galaxy screenshot of heatmapper usage and output" />
</a>

___________________________________
<a name="multiprofiler"></a>
#### How can I compare the average signal for X- and autosomal genes for 2 or more different sequencing experiments?

Make sure you're familiar with computeMatrix and profiler before using this protocol.

* tools:
    * Filter data on any column using simple expressions
    * computeMatrix
    * profiler
    * (plotting the summary plots for multiple samples)

* inputs:
    * several bigWig files (one for each sequencing experiment you would like to compare)
    * two BED files, one with X-chromosomal and one with autosomal genes
___________________________________

##### How to obtain a BED file for X chromosomal and autosomal genes each

1. download a full list of genes via "Get Data" &rarr; "UCSC main table browser" &rarr; group:"Genes and Gene Predictions" &rarr; tracks: (e.g.) "RefSeqGenes" &rarr; send to Galaxy

2. filter the list twice using the tool __"Filter data on any column using simple expressions"__ 

    - first use the expression: c1=="chrX" to filter the list of all genes &rarr; this will generate a list of X-linked genes
    - then re-run the filtering, now with c1!="chrX", which will generate a list of genes that do not belong to chromosome X (!= indicates "not matching")

##### Compute the average values for X and autosomal genes 

* use [computeMatrix][] for all of the signal files (bigWig format) at once

    * supply both filtered BED files (click on "Add new regions to plot" once) and label them
    * indicate the corresponding signal files

* now use [profiler][] on the resulting file

    * important: display the "advanced output options" and select "save the data underlying the average profile" &rarr; this will generate a table in addition to the summary plot images

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_profiles_XvsA02.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_profiles_XvsA02.png" Title="combined profiles for different histone marks" />
</a>

------------------------------------
[multiBamCoverage]: https://github.com/fidelram/deepTools/wiki/QC
[plotFingerprint]: https://github.com/fidelram/deepTools/wiki/QC
[computeGCBias]: https://github.com/fidelram/deepTools/wiki/QC
[bamCoverage]: https://github.com/fidelram/deepTools/wiki/Normalizations
[bamCompare]: https://github.com/fidelram/deepTools/wiki/Normalizations
[correctGCbias]: https://github.com/fidelram/deepTools/wiki/Normalizations
[computeMatrix]: https://github.com/fidelram/deepTools/wiki/Visualizations
[heatmapper]: https://github.com/fidelram/deepTools/wiki/Visualizations
[profiler]: https://github.com/fidelram/deepTools/wiki/Visualizations 

[Galaxy]: http://galaxyproject.org/ "General Galaxy platform from Penn State"
[GEO]: http://www.ncbi.nlm.nih.gov/geo/ "GEO database"
[Roadmap project]: http://www.roadmapepigenomics.org/data "Roadmap web site"
[UCSC]: http://genome.ucsc.edu/ "UCSC Genome web site"
[BioMart]: http://www.biomart.org/ "Biomart web site"
[deepTools Galaxy]: http://deeptools.ie-freiburg.mpg.de/ "deepTools Galaxy at the Max-Planck-Institute of Immunobiology and Epigenetics"

[2bit]: https://github.com/fidelram/deepTools/wiki/Glossary#wiki-2bit "binary file for storage of genome sequences"
[BAM]: https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bam "binary version of a SAM file; contains all information about aligned reads"
[bed]: https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bed "text file that usually contains gene information such as chromosome, gene start, gene end, gene name, strand information - can be used for any genomic region representation"
[BED]: https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bed "text file that usually contains gene information such as chromosome, gene start, gene end, gene name, strand information - can be used for any genomic region representation"
[bedGraph]: https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bedgraph "text file that contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[bigWig]: https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bigwig "binary version of a bedGraph file; contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[FASTA]: https://github.com/fidelram/deepTools/wiki/Glossary#wiki-fasta "simple text-file containing nucleotide or protein sequences"
[FASTQ]: https://github.com/fidelram/deepTools/wiki/Glossary#wiki-fastq "text file of raw reads (almost straight out of the sequencer)"
[SAM]: https://github.com/fidelram/deepTools/wiki/Glossary#wiki-sam "text file containing all information about aligned reads"
[bin]: https://github.com/fidelram/deepTools/wiki/Glossary#terminology "typically a small region of the genome, used to 'store' a score; created by artificially dividing the genome"
[read]: https://github.com/fidelram/deepTools/wiki/Glossary#terminology "the DNA piece that was actually sequenced  ("read") by the sequencing machine (usually between 30 to 100 bp long, depending on the read-length of the sequencing protocol)" 
[input]: https://github.com/fidelram/deepTools/wiki/Glossary#terminology "confusing, albeit commonly used name for the 'no-antibody' control sample for ChIP experiments"

