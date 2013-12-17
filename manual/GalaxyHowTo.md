How do I...?
======================

This site is meant to give you concise guidance to specific tasks you may want to perform. There are many more ways in which you can use [deepTools Galaxy][] than those described here, so be creative once you're comfortable with using them.

* [I have downloaded/received a FASTQ file - how do I generate a file I can look at in a Genome Browser?](#FASTQ2IGV)
* [How can I assess the reproducibility of my sequencing replicates?](#repCorr)
* [How do I get an input-normalized ChIP-seq coverage file?](#InputNorm)
* [How do I know whether my sample is GC biased? And if yes, how do I correct for it?](#GC)
* [How do I get a (clustered) heatmap of sequencing-depth-normalized read coverages around the transcription start site of all genes?](#HM)
* [How can I compare the average signal for X- and autosomal genes for 2 or more different sequencing experiments](#profiler)

[Back to the general Galaxy help page](https://github.com/fidelram/deepTools/blob/master/manual/GalaxyHelp.md)
------------------------------------------------------------------------------------------

#### I have downloaded/received a [FASTQ][] file - how do I generate a file I can look at in a Genome Browser?<a name="FASTQ2IGV">
#### How can I assess the reproducibility of my sequencing replicates?<a name="repCorr">
#### How do I get an input-normalized ChIP-seq coverage file?]<a name="InputNorm">
#### How do I know whether my sample is GC biased? And if yes, how do I correct for it?<a name="#GC")
#### How do I get a (clustered) heatmap of sequencing-depth-normalized read coverages around the transcription start site of all genes?<a name="#HM">
#### How can I compare the average signal for X- and autosomal genes for 2 or more different sequencing experiments](#profiler)


-----------------------------------------------------------------------------------
[Back to the general Galaxy help page](https://github.com/fidelram/deepTools/blob/master/manual/GalaxyHelp.md)
-----------------------------------------------------------------------------------
[Galaxy]: http://galaxyproject.org/ "General Galaxy platform from Penn State"
[GEO]: http://www.ncbi.nlm.nih.gov/geo/ "GEO database"
[Roadmap project]: http://www.roadmapepigenomics.org/data "Roadmap web site"
[UCSC]: http://genome.ucsc.edu/ "UCSC Genome web site"
[BioMart]: http://www.biomart.org/ "Biomart web site"
[deepTools Galaxy]: http://deeptools.ie-freiburg.mpg.de/ "deepTools Galaxy at the Max-Planck-Institute of Immunobiology and Epigenetics"

[MACS]: http://www.ncbi.nlm.nih.gov/pubmed/22936215 "How to use MACS, Nature Protocols"
[CCAT]: http://www.ncbi.nlm.nih.gov/pubmed/20371496 "CCAT original publication"
[SICER]: http://bioinformatics.oxfordjournals.org/content/25/15/1952.full "SICER original publication"

[BAM]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "binary version of a SAM file; contains all information about aligned reads"
[SAM]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file containing all information about aligned reads"
[bigWig]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "binary version of a bedGraph file; contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[bedGraph]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file that contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[FASTQ]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file of raw reads (almost straight out of the sequencer)"
