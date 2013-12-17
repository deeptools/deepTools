###[Tipps](#Tipps)

* [How to optimize the computation time](#compTime)
* [How to integrate the deepTools results with other downstream analyses (outside of Galaxy)](#integrate)

###[FAQ](#FAQ)
* [I've reached my quota - what can I do to save some space?](#quota)
* [How can I determine basic parameters of a BAM file?](#BAMparams)
* [When should I exclude regions from computeGCbias?](#excludeGC)
* [How can I increase the resolution of the heatmap?](#hmresolution)


----------------------------------------------------------------------------------------------------
<a name="Tipps">
Tipps for optimal deepTools usage
==========================================


#### Optimizing the computation time <a name="compTime">
* when you're playing around with the tools to see what kinds of results they will produce: limit the operation to one chromosome only to __save computation time__! ("advanced output options" --> "Region of the genome to limit the operation to")

#### Integrate the deepTools results with other downstream analyses (outside of Galaxy) <a name="integrate">
* you can save all the data tables underlying every image produced by deepTools, i.e. if you would like to plot the average profiles in a different way, you could download the corresponding data (after ticking the profiler option at "advanced output options") and import them into R, Excel, GraphPadPrism etc.


<a name="FAQ">
Frequently asked questions
===========================

#### I've reached my quota - what can I do to save some space? <a name="quota">
1. make sure that all the data sets you deleted are __permanently__ eliminated from our disks: go to the history option button and select "Purge deleted data sets", then hit the "refresh" button on top of your history panel
2. download all data sets for which you've completed the analysis, then remove the data sets (click on the "x" and then make sure they're purged (see above)


#### How can I determine basic parameters of a BAM file, such as the number of reads, read length, duplication rate and average DNA fragment length? <a name="BAMparams">

Eventhough [MACS][] is meant to do peak calling for you, it also outputs a number of useful information such as those listed above.
Simply run MACS on the BAM file that you would like to gain the information for and check the .xls file from the MACS output. It will list:

* tag length = read length
* duplication rate
* number of tags = number of reads
* d = distance = average DNA fragment size


#### When should I exclude regions from computeGCbias? <a name="excludeGC">

In general, we recommend that you should only correct for GC bias (using computeGCbias followed by correctGCbias) if you observe that the majority of the genome (the region between 30-60%) is continuously GC-biased __and__ you want to compare this sample with another sample that is not GC-biased.

Sometimes, a certain GC bias is expected, for example for ChIP samples of H3K4me3 in mammalian samples where GC-rich promoters are expected to be enriched. To not confound the GC bias caused by the library preparation with the inherent, expected GC bias, we incorporated the possibility to supply a file of regions to computeGCbias that will be excluded from the GC bias calculation. This file should typically contain those regions that one expects to be significantly enriched per se. This way, the computeGCbias will focus on background regions.


#### The heatmap I generated looks very "coarse", I would like a much more fine-grained image. <a name="hmresolution">
* decrease the __bin size__ when generating the matrix using computeMatrix
  * go to "advanced options" --> "Length, in base pairs, of the non-overlapping bin for averaging the score over the regions length" --> define a smaller value, e.g. 50 or 25 bp

* make sure, however, that you used a sufficiently small bin size when calculating the bigWig file, though (if generated with deepTools, you can check the option "bin size")
 
