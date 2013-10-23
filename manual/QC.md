QC of aligned reads
====================

## Table of Contents  
  * [bamCorrelate](#bamCorrelate)
  * [computeGCbias](#computeGCbias)
  * [bamFingerprint](#bamFingerprint)


<a name="bamCorrelate"/>
bamCorrelate
-------------


This tool is useful to assess the overall similarity of different BAM
files. A typical application is to check the correlation between
replicates or published data sets.

### What it does
The tool splits the genomes into bins of a given length. For each bin,
the number of reads found in each BAM file is counted and a
correlation is computed for all pairs of BAM files.

### Important parameters
bamCorrelate can be run in 2 modes: _bins_ and _bed_.

In the bins mode, the correlation is computed based on equal length bins. The user has to specifcy the _number of bins_. This is useful to assess the overall similarity of BAM files,
but due to the random sampling of the genome, outliers might be produced that will skew the correlation values.

In the BED-file options, the user has to supply a list of genomic regions in BED format (chr, start, end) in addition to (a) BAM file(s). bamCorrelate subsequently only calculates the
number of overlapping reads for these regions. This can be used, for example, to compare the ChIP-seq coverages of two different samples in peak regions only.

### Output files:
    - plot
    - matrix

### Example Figures

replicates, input vs ChIP, patient 1 vs patient2



<a name="computeGCbias"/>
computeGCbias
--------------

This tool computes the GC bias using the method proposed by [Benjamini
and Speed][] (see below for more explanations). 


### What it does (reference to Benjamini)

The basic assumption of the GC bias diagnosis is that the distribution of sequenced reads might not be
uniform across the genome, i.e. not all regions of the genome are similarly well sequenced.

computeGCbias estimates how many reads with what kind of GC content one
should have sequenced. This is based on the methods published by [Benjamini and Speed][].
The tool first determines how many regions the specific reference genome contains for each amount of GC content,
i.e. how many regions in the genome have 50% GC (or 10% GC or 90%
GC or...).  For this, it samples a large number of equally sized genome bins
and counts how many times we see a bin with 50% GC (or 10% GC or 90%
or...). These __expected values are independent of any sequencing, but they do depend on the respective reference genome__.
This means, that the expected values will of course differ between mouse and fruit fly due to their genome's
different GC contents, but it also means that strong biases in the reference genome assembly might lead to a false positive diagnosis of GC bias.

After the expected values, the tool samples the BAM file of sequenced reads. Instead of noting how many genomic regions
there are per GC content, we now count the __reads per GC content__. 

### Output files

  + __Diagnostic plot__
    - box plot of absolute read numbers per genomic GC bin
    - x-y plot of observed/expected read ratios per genomic GC content bin
    
  + __Data matrix__
    - to be used for GC correction with _correctGCbias_
    
### What the plots tell you
In an ideal sample without GC bias, the ratio of observed/expected values
should be close to 1 for all GC content bins.

However, due to PCR (over)amplifications, the majority of ChIP samples usually shows a
significant bias towards reads with high GC content (>50%) and a depletion of reads from GC-poor regions.

### Example figures
  + Figure: strong bias vs. no bias vs. histone mark with inherent GC bias, GC bias of simulated data



<a name="bamFingerprint"/>
bamFingerprint
---------------

### What it does
This tool is based on a method developed by [Diaz et al.][].
For factors that will enrich well-defined, rather narrow regions (e.g. transcription factors such as p300), the resulting plot can be used to assess
the strength of a ChIP, i.e. whether the signal of the enrichment can be clearly distinguished from the background.

The tool first samples indexed [BAM][] files and counts all reads overlapping
a window (bin) of specified length. These counts are then sorted
according to their rank and the cumulative sum of read counts is
plotted.

###	Output files: 
  + plot
  + matrix

### What the plots tell you
An ideal input with perfect uniform distribution of reads
along the genome (i.e. without enrichments in open chromatin etc.)
should generate a straight diagonal line. A very specific and strong
ChIP enrichment will be indicated by a prominent and steep rise of the
cumulative sum towards the highest rank. This means that a big chunk
of reads from the ChIP sample is located in few bins which corresponds
to high, narrow enrichments seen for transcription factors.

  + Figure: good TF vs. bad TF vs. broad histone mark


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
