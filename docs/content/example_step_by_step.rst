Step-by-step protocols
========================

.. contents:: 
    :local:


How can I do...? <a name="HowTo"></a>
--------------------------------------------------

This section is meant to give you quick guidance to specific tasks you may want to perform. We're using screenshots from Galaxy here, if you're using the command-line version, you can easily follow the given examples by typing the program name and the help option (e.g. /deepTools/bin/bamCoverage --help) which will show you all the parameters and options, most of them named very similarly to those in Galaxy.

For each "recipe" here, you will find the screenshot of the tool and the input parameters on the left hand side (we marked non-default, _user-specified entries_) and screenshots of the output on the right hand side. Do let us know if you spot things that are missing, should be explained better or are plain confusing!

There are many more ways in which you can use [deepTools Galaxy][] than those described here, so be creative once you're comfortable with using them. For detailed explanations of what the tools do, follow the links.

> All recipes assume that you have uploaded your files into a Galaxy instance with a deepTools installation, e.g. [deepTools Galaxy][].

> If you would like to try out the protocols with __sample data__, go to [deepTools Galaxy][]  &rarr; "Shared Data"  &rarr; "Data Libraries"  &rarr; "Sample Data". Use one file from each folder. E.g., import 1 BAM file from the folder "mapped reads", 1 bigwig file from the folder "normalized read coverages" (preferably "Dmel_log2ratio....bigwig", and the .bed file "Drosophila housekeeping genes" from the folder "annotation data" into your current Galaxy history. For testing our protocols via the command line, you can download the sample files to your computer by clicking on the triangle right next to the file name.

___________________________________
<a name="FASTQ2IGV"></a>
#### I have downloaded/received a [BAM][] file - how do I generate a file I can look at in a Genome Browser?

* tool: [bamCoverage][]
* input: your BAM file

Note: BAM files can also be viewed in Genome Browsers, however, they're large and tend to freeze the applications. Generating bigWig files of read coverages will help you a lot in this regard. In addition, if you have more than one sample you'd like to look at, it is helpful to normalize all of them to 1x sequencing depth.

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamCoverage.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamCoverage.png" Title="deepTools Galaxy screenshot of bamCoverage usage and output" />
</a>

___________________________________
<a name="repCorr"></a>
#### How can I assess the reproducibility of my sequencing replicates?

* tool: [bamCorrelate][]
* input: BAM files
    * you can compare as many samples as you want - the more you put at the same time, the longer the computation takes

* output: heatmap of correlations - the closer two samples are to each other, the more similar their read coverages

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamCorrelate.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamCorrelate.png" Title="deepTools Galaxy screenshot of bamCorrelate usage and output" />
</a>

___________________________________
<a name="GC"></a>
#### How do I know whether my sample is GC biased? And if yes, how do I correct for it?

* you need a BAM file of your sample in question
* use the tool [computeGCbias][] on that BAM file (default settings, just make sure your reference genome and genome size are matching)

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_computeGCbias.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_computeGCbias.png" Title="deepTools Galaxy screenshot of computeGCbias" />
</a>


* have a look at the image that is produced and compare it to the examples [here](https://github.com/fidelram/deepTools/wiki/QC)
* if your sample shows an almost linear increase in exp/obs coverage (on the log scale of the lower plot), then you should consider correcting the GC bias - _if_ you think that the biological interpretation of this data would otherwise be compromised (e.g. by comparing it to another sample that does not have an inherent GC bias)

    + the GC bias can be corrected with the tool [correctGCbias][] using the second output of the computeGCbias tool that you had to run anyway

    + CAUTION!! correctGCbias will add reads to otherwise depleted regions (typically GC-poor regions), that means that you should __not__ remove duplicates in any downstream analyses based on the GC-corrected BAM file (we therefore recommend to remove duplicates before doing the correction so that only those duplicate reads are kept that were produced by the GC correction procedure)

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_correctGCbias.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_correctGCbias.png" Title="deepTools Galaxy screenshot of correctGCbias usage and output" />
</a>

___________________________________
<a name="InputNorm"></a>
#### How do I get an input-normalized ChIP-seq coverage file?

* input: you need two BAM files, one for the input, one for the ChIP-seq experiment
* tool: [bamCompare][] with ChIP = treatment, input = control sample

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamCompare.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamCompare.png" Title="deepTools Galaxy screenshot of bamCompare usage and output" />
</a>

___________________________________
<a name="fprint"></a>
#### How can I compare the ChIP strength for different ChIP experiments?

* tool: [bamFingerprint][]
* input: as many BAM files as you'd like to compare. Make sure you get all the labels right!

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamFingerprint.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_bamFingerprint.png" Title="deepTools Galaxy screenshot of bamFingerprint" />
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
    * if you would like to cluster the signals, choose "kmeans clustering" (last option of "advanced options") with a reasonable number of clusters (usually between 2 to 7)

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
    - then re-run the filtering, now with c1!="chrX" which will generate a list of genes that do not belong to chromosome X (!= indicates "not matching")

##### Compute the average values for X and autosomal genes 

* use [computeMatrix][] for __each__ signal file (bigWig) (you only need to specify all the parameters once, then use the re-run button underneath the first data set and just replace the signal file with the next one)

    * supply both filtered BED files (click on "Add new regions to plot" once) and label them
    * indicate the corresponding signal file
    * make sure to __re-name__ every data set in the history once computeMatrix is done so that you can easily keep track of which matrix was based on which bigWig file (you can always find these information by clicking on the i-button in the respective data set)

* now use [profiler][] for every file you generated with computeMatrix

    * important: display the "advanced output options" and select "save the data underlying the average profile" &rarr; this will generate a table in addition to the summary plot images

* now you have at least 2 separate images of profiles - one for each bigWig file - you can either leave it like this or use another script that will plot all the summary plots in one image at once

    * this tool is called "Plotting the summary plots for multiple signals"
    * it uses the tables generated by profiler
    * for each group of genes (in this case, X and autosomal genes = 2 groups), you can assign a color

The result could look like this:

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_profiles_XvsA02.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_profiles_XvsA02.png" Title="combined profiles for different histone marks" />
</a>

As you may have noticed, this task requires several steps that are repeated. Here is a screenshot of how the Galaxy workflow would look like (you can find it under "Shared Data" &rarr; "Published Workflows" &rarr; "Summary plots for X and autosomal genes" where we have constructed it with the example histone marks from the Data Library. Be aware that running this workflow will take up quite some computation timing, but it won't require much input from your part - so start if before you go off for lunch...)

If you're not sure how to use the published workflow, please read [this entry](#workflow) or go to the central [Galaxy learning page full of tutorials]( https://wiki.galaxyproject.org/Learn/Screencasts#Tutorials "Galaxy Tutorials").

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_profiles_XvsA.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_profiles_XvsA.png" Title="Screenshot of the workflow designed for the above described task" />
</a>

------------------------------------
[bamCorrelate]: https://github.com/fidelram/deepTools/wiki/QC
[bamFingerprint]: https://github.com/fidelram/deepTools/wiki/QC
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

