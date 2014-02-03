This document contains some chapters of our wiki on deepTools usage for NGS data analysis. For the most updated version of our help site and for **more information about deepTools, a brief introduction into Galaxy** as well as **step-by-step protocols**, please visit: https://github.com/fidelram/deepTools/wiki

We included many screenshots of our Galaxy deepTools web server to explain the usage of our tools. In addition, we show the commands for the stand-alone usage, as they often indicate the options that one should pay attention to more succinctly.  

Why we built deepTools <a name="why"></a>
-------------------------------------------

The main reason why deepTools was started is the simple fact that in 2011 we could not find tools that met all our needs for NGS data analysis. While there were individual tools for separate tasks, we wanted software that would fulfill *all* of the following criteria:

* **efficiently extract reads from BAM files** and perform various computations on them
* **turn BAM files of aligned reads into bigWig files** using different normalization strategies
* make use of **multiple processors** (speed!)
* generation of **highly customizable images** (change colors, size, labels, file format etc.)
* enable **customized down-stream analyses** which requires that every data set that is being produced can be stored by the user
* **modular approach** - compatibility, flexibility, scalability (i.e. we can add more and more modules making use of established methods)


The flow chart below depicts the different tool modules that are currently available within deepTools (deepTools modules are written in bold red and black font).

![flowChartI](https://raw.github.com/fidelram/deepTools/master/examples/flowChart_BAMtoBIGWIG_small.png "Average analysis and QC workflow")



[Benjamini and Speed]: http://nar.oxfordjournals.org/content/40/10/e72 "Nucleic Acids Research (2012)"
[Diaz et al.]: http://www.degruyter.com/view/j/sagmb.2012.11.issue-3/1544-6115.1750/1544-6115.1750.xml "Stat. Appl. Gen. Mol. Biol. (2012)"
[IGV]: http://www.broadinstitute.org/igv/ "Integrative Genome Browser developed by the Broad Institute"
[bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml


 
