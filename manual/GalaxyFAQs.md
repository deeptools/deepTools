Frequently asked questions
===========================

#### [Galaxy-specific questions](#GalSpecific)
* [I've reached my quota - what can I do to save some space?](#quota)
* [What is the best way to integrate the deepTools results with other downstream analyses (outside of Galaxy)?](#integrate)
* [How can I determine basic parameters of a BAM file?](#BAMparams)

#### [General deepTools-related questions](#general)
* [I just want to try out a tool, how can I optimize the computation time?](#compTime)
* [When should I exclude regions from computeGCbias?](#excludeGC)
* [Does it speed up the computation if I limit bamCorrelate to one chromosome, but keep the same numbers and sizes of sampling bins?](#bamCorrelateLimit)
* [How do I calculate the effective genome size for an organism that's not in your list?](#effGenomeSize)
* [How can I increase the resolution of the heatmap?](#hmresolution)


----------------------------------------------------------------------------------------------------

Galaxy-specific questions <a name="GalSpecific">
--------------------------------------------------

#### I've reached my quota - what can I do to save some space? <a name="quota">
1. make sure that all the data sets you deleted are __permanently__ eliminated from our disks: go to the history option button and select "Purge deleted data sets", then hit the "refresh" button on top of your history panel
2. download all data sets for which you've completed the analysis, then remove the data sets (click on the "x" and then make sure they're purged (see above)


#### What's the best way to integrate the deepTools results with other downstream analyses (outside of Galaxy) <a name="integrate">
* you can save all the data tables underlying every image produced by deepTools, i.e. if you would like to plot the average profiles in a different way, you could download the corresponding data (after ticking the profiler option at "advanced output options") and import them into R, Excel, GraphPadPrism etc.


#### How can I determine basic parameters of a BAM file, such as the number of reads, read length, duplication rate and average DNA fragment length? <a name="BAMparams">

Eventhough [MACS][] is meant to do peak calling for you, it also outputs a number of useful information such as those listed above.
Simply run MACS on the BAM file that you would like to gain the information for and check the .xls file from the MACS output. It will list:

* tag length = read length
* duplication rate
* number of tags = number of reads
* d = distance = average DNA fragment size



General deepTools-related questions <a name="general">
--------------------------------------------------------------

#### How can I test a tool with little computation time? <a name="compTime">
* when you're playing around with the tools to see what kinds of results they will produce: limit the operation to one chromosome only to __save computation time__! ("advanced output options" --> "Region of the genome to limit the operation to")


#### When should I exclude regions from computeGCbias? <a name="excludeGC">

In general, we recommend that you should only correct for GC bias (using computeGCbias followed by correctGCbias) if you observe that the majority of the genome (the region between 30-60%) is continuously GC-biased __and__ you want to compare this sample with another sample that is not GC-biased.

Sometimes, a certain GC bias is expected, for example for ChIP samples of H3K4me3 in mammalian samples where GC-rich promoters are expected to be enriched. To not confound the GC bias caused by the library preparation with the inherent, expected GC bias, we incorporated the possibility to supply a file of regions to computeGCbias that will be excluded from the GC bias calculation. This file should typically contain those regions that one expects to be significantly enriched per se. This way, the computeGCbias will focus on background regions.


#### Does it speed up the computation if I limit bamCorrelate to one chromosome, but keep the same numbers and sizes of sampling bins?<a name="bamCorrelateLimit">
Yes. However, the way bamCorrelate (and all the other deepTools handle the option "limit the computation to a specific region" is as follows: first, the _entire_ genome represented in the BAM file will be regarded and sampled, _then_ all the regions or sampled bins that do not overlap with the region indicated by the user will be discarded. This means that if you wanted 10,000 bins to be sampled and you focus on, let's say, chromosome 2, the final computation will not be performed on the whole set of 10,000 bins, but only on those bins that overlap with chromosome 2.

#### How do I calculate the effective genome size for an organism that's not in your list?]<a name="effGenomeSize">
This is something you will have to find a solution outside of deepTools at the moment. We suggest to run faCount from UCSC tools. If you used multi-read alignment (e.g. with bowtie2), then you can use that tool to report the total number of bases as well as the number of unmapped bp, indicated by 'N'. The effective genome size is the total number of reads minus the number of 'N'.

#### The heatmap I generated looks very "coarse", I would like a much more fine-grained image. <a name="hmresolution">
* decrease the __bin size__ when generating the matrix using computeMatrix
  * go to "advanced options" --> "Length, in base pairs, of the non-overlapping bin for averaging the score over the regions length" --> define a smaller value, e.g. 50 or 25 bp

* make sure, however, that you used a sufficiently small bin size when calculating the bigWig file, though (if generated with deepTools, you can check the option "bin size")
 
