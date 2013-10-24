Normalization of BAM files
===========================

deepTools contains 3 tools for the normalization of [BAM][] files:

1. _correctGCbias_: in case, you would like to normalize your read distributions to fit the expected GC values, you can use the output from [computeGCbias](https://raw.github.com/fidelram/deepTools/master/manual/QC.md "go to the chapter about data QC") and produce a GC-corrected [BAM]-file.
2. _bamCoverage_: this tool converts a single [BAM][] file into a [bigWig][] file, enabling you to normalize for sequencing depth.
3. _bamCompare_: like bamCoverage, this tool produces a normalized [bigWig][] file, but it takes 2 [BAM][] files, normalizes them for sequencing depth and subsequently performs a mathematical operation of your choice, i.e. it can output the ratio of the read coverages in both files or the like.


[Here](https://docs.google.com/file/d/0B8DPnFM4SLr2UjdYNkQ0dElEMm8/edit?usp=sharing "How to get from aligned reads to coverage profiles using deepTools") you can download slides that we used for teaching. They contain additional details about how the coverage files are generated and normalized.


![bamToBigWig](https://raw.github.com/fidelram/deepTools/master/examples/norm_IGVsnapshot_indFiles.png "snapshots of bigWig files loaded in IGV")

## Table of Content

  * [correctGCbias](#correctGCbias)
  * [bamCoverage](#bamCoverage)
  * [bamCompare](#bamCompare)



<a name="correctGCbias"/>
CorrectGCbias
---------------

### What it does (uses output from computeGCbias)
This tool requires the output from computeGCBias to correct the given
[BAM][] files according to the method proposed by [Benjamini and Speed][].  The resulting [BAM][] files can be used in any
downstream analyses, but __be aware that you should not filter out
duplicates from here on__.

### output
 + __GC-normalized BAM file__


<a name="bamCoverage"/>
bamCoverage
------------

### What it does
Given a BAM file, this tool generates a [bigWig][] or [bedGraph][] file of fragment or read coverages. The way the method works is by first calculating all the number of reads (either extended to match the fragment length or not) that overlap each bin in the genome. Bins with zero counts are skipped, i.e. not added to the output file. The resulting read counts can be normalized using either a given scaling factor, the RPKM formula or to get a 1x depth of coverage (RPGC).

### output
  + __coverage file__ either in [bigWig][] or [bedGraph][] format


### Usage

Here's an exemplary command to generate a single [bigWig][] file out of a single [BAM][] file:

    $/deepTools-1.5/bin/bamCoverage -b corrected_counts.bam -bs 10 --normalizeTo1x 2150570000 --fragmentLength 200 -o Coverage.GCcorrected.SeqDepthNorm.bw --ignoreForNormalization chrX

  + The bin size __(-bs)__ can be chosen completely to your liking. The smaller it is, the bigger your file will be.
  + This was a mouse sample, therefore the effective genome size for mouse had to be indicated once it was decided that the file should be normalize to 1x coverage.
  + Chromosome X was excluded from sampling the regions for normalization as the sample was from a male mouse that therefore contained pairs of autosome, but only a single X chromosome.
  + The fragment length of 200 bp is only the fall-back option of bamCoverage as the sample provided here was done with paired-end sequencing. Only in case of singletons will bamCoverage resort to the user-specified fragment length.

##### Some (more) parameters to pay special attention to

  + --ignoreDuplicates - important! in case where you normalized for GC bias using correctGCbias, you should absolutely  __NOT__ set this parameter
 

<a name="bamCompare"/>
bamCompare
------------

### What it does

This tool compares two BAM files based on the number of mapped
reads. To compare the BAM files, the genome is partitioned into bins
of equal size, the reads are counted for each bin and each BAM file
and finally, a summarizing value is reported.  This value can be the
ratio of the number of reads per bin, the log2 of the ratio or the
difference.  This tool can normalize the number of reads on each BAM
file using the SES method proposed by [Diaz et al.][] Normalization based on read counts is also
available. If paired-end
reads are present, the fragment length reported in the BAM file is
used by default.

### output file
  + same as for bamCoverage, except that you now obtain __1__ coverage file that is based on __2__ [BAM][] files.

### Usage

Here's an example command that generated the log2(ChIP/Input) values.

    $ /deepTools-1.5/bin/bamCompare -b1 ChIP.bam -b2 Input.bam -bs 25 -f 200 --missingDataAsZero no --ratio log2 --scaleFactorsMethod SES -o log2ratio_ChIP_vs_Input.bw
    
  + like for bamCoverage, the bin size is completely up to the user
  + the fragment size (-f) will only be taken into consideration for reads without mates
  + the SES method was used for normalization as the ChIP sample was done for a histone mark with highly localized enrichments (similar to the left-most plot of the [fingerprint-examples](https://github.com/fidelram/deepTools/blob/master/manual/QC.md "change to QC chapter to have a look at the plots")

##### Some (more) parameters to pay special attention to

 + --scaleFactorsMethod - here you can choose how you would like to normalize to account for variation in sequencing depths. We provide the simple normalization total read count or the more sophisticated signal extraction (SES) method proposed by [Diaz et al.][]. We recommend to use SES only for those cases where yhe distinction between input and ChIP is very clear in the [bamFingerprint plots](https://github.com/fidelram/deepTools/blob/master/manual/QC.md "go back to the QC chapter"). This is usually the case for transcription factors and sharply defined histone marks such as H3K4me3.
  + --ratio - here you get to choose how you want the two input files to be compared, e.g. by taking the ratio or by subtracting b2 from b1 etc. In case you do want to subtract one sample from the other, you will have to choose whether you want to normalize to 1x coverage (--normalizeTo1x) or to __r__eads __p__er __k__ilobase (--normalizeUsingRPKM; similar to RNA-seq normalization schemes)
 + --ignoreForNormalization - same as for bamCoverage


![bamCompare](https://raw.github.com/fidelram/deepTools/master/examples/norm_bamCompare.png "Mathematical operations for comparing 2 BAM files implemented in bamCompare")



-----------------------------------------------------------------------------------
[BAM]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "binary version of a SAM file; contains all information about aligned reads"
[SAM]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file containing all information about aligned reads"
[bigWig]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "binary version of a bedGraph file; contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[bedGraph]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file that contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[FASTQ]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file of raw reads (almost straight out of the sequencer)"
### References
[Benjamini and Speed]: http://nar.oxfordjournals.org/content/40/10/e72 "Nucleic Acids Research (2012)"
[Diaz et al.]: http://www.degruyter.com/view/j/sagmb.2012.11.issue-3/1544-6115.1750/1544-6115.1750.xml "Stat. Appl. Gen. Mol. Biol. (2012)"


This tool is developed by the [Bioinformatics Facility](http://www1.ie-freiburg.mpg.de/bioinformaticsfac) at the [Max Planck Institute for Immunobiology and Epigenetics, Freiburg](http://www1.ie-freiburg.mpg.de/).
