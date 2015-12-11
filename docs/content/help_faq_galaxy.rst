Galaxy-related FAQ
===================

.. contents:: 
    :local:

I've reached my quota - what can I do to save some space?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1. make sure that all the data sets you deleted are **permanently** eliminated from our disks: go to the history option button and select "Purge deleted data sets", then hit the "refresh" button on top of your history panel
2. download all data sets for which you've completed the analysis, then remove the data sets (click on the "x" and then **make sure they're purged** (see above))

Copying from one history to another doesn't work for me - the data set simply doesn't show up in the target history!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you've copied a data set from one history to another, check two things:
* do you see the destination history in your history panel, i.e. does the title of the current history panel match the name of the destination history you selected in the main frame?
* hit the refresh button

<a href="https://raw.github.com/fidelram/deepTools/master/examples/Gal_historyReload.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/Gal_historyReload.png" Title="Galaxy history refresh button" />
</a>

How can I use a published workflow?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You **must register** if you want to use the workflows within [deepTools Galaxy][]. ("User" &rarr; "Register" - all you have to supply is an email address)

You can find workflows that are public or specifically shared with you by another user via "Shared Data" &rarr; "Published Workflows". Click on the triangle next to the workflow you're interested in and select "import".

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_wf01.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_wf01.png" Title="Finding published workflows" />
</a>

A green box should appear, there you select "start using this workflow" which should lead you to your own workflow menu (that you can always access via the top menu "Workflow"). Here, you should now see a workflow labeled "imported: ....". If you want to use the workflow right away, click on the triangle and select "Run". The workflow should now be available within the Galaxy main data frame and should be waiting for your input.

<a href="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_wf02.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/GalHow_wf02.png" Title="Finding published workflows" />
</a>

I would like to use one of your workflows - not in the deepTools Galaxy, but in the local Galaxy instance provided by my institute. Is that possible?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yes, it is possible. The only requirement is that your local Galaxy has a recent installation of deepTools.

Go to the workflows, click on the ones you're interested in and go to "Download". This will save the workflows into .ga files on your computer. Now go to your local Galaxy installation and login. Go to the workflow menu and select "import workflow" (top right hand corner of the page). Click on "Browse" and select the saved workflow. If you have the same tool versions installed in your local Galaxy, these workflows should work right away.

How can I have a look at the continuous read coverages from bigWig files? Which Genome Browser do you recommend?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are 2 popular Genome Browsers to visualize continuous data: [UCSC](http://genome.ucsc.edu/cgi-bin/hgGateway?redirect=manual&source=genome-euro.ucsc.edu "UCSC Genome Browser") and [IGV](http://www.broadinstitute.org/igv/ "IGV Genome Browser").

IGV (recommended)
""""""""""""""""""

We recommend to download the [IGV Genome Browser](http://www.broadinstitute.org/igv/ "IGV Genome Browser"). Any user working for academia will be granted access after an informal registration. IGV itself needs an up-to-date Java installation and considerable amounts of RAM. It's usage is rather intuitive and the display can be easily customized. In addition, you can download genome-wide annotation data that can be displayed together with your own data.

To display data in IGV, do the following steps:

1. Go to http://www.broadinstitute.org/igv/, register and download IGV
2. Unpack the IGV archive and change to the extracted IGV folder
3. Use the igv.bat (Windows), igv.sh (Linux) or igv.command (OSX) to start IGV (for more information please read the included readme.txt file or the IGV documentation)
4. Choose the genome version of the file(s) you would like to visualize (e.g. dm3) THIS IS THE MOST IMPORTANT STEP! IGV will not detect the genome version automatically, i.e. if you select mm9 but your file is based on human data, it will still be displayed without an error message (but with the wrong positions, obviously!)
5. Go to your deepTools Galaxy server (http://deeptools.ie-freiburg.mpg.de/) and navigate to your data set of choice
6. Click on your data set so that you see its details like in  the screenshot below (_Keep in mind that not all datasets can be visualized in IGV or UCSC._ We recommend to use [bigWig][] or [BED][] files for visualization.)


<a href="https://raw.github.com/fidelram/deepTools/master/examples/Gal_FAQ_IGV_dataset.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/Gal_FAQ_IGV_dataset.png" Title="Screenshot of a bigWig file in Galaxy" />
</a>

Now click on “display with IGV local” to visualize your data set in IGV that should already be running on your computer. _(“display with IGV Web current” can be used if you do not have an installed IGV. It will start an IGV web start version. We do not recommend that option.)_

Here's a screenshot of a typical bigWig file display:
<a href="https://raw.github.com/fidelram/deepTools/master/examples/Gal_FAQ_IGV.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/Gal_FAQ_IGV.png" Title="Screenshot of IGV browser display of bigWig files" />
</a>


For more information, check out the [IGV documentation](http://www.broadinstitute.org/software/igv/UserGuide "IGV User Guide").

UCSC
""""""

There is a direct link from within deepTools Galaxy to stream a data set to UCSC. You can find it in the data set tiles: "display at UCSC", like here:

<a href="https://raw.github.com/fidelram/deepTools/master/examples/Gal_FAQ_UCSC_dataset.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/Gal_FAQ_UCSC_dataset.png" Title="Screenshot of bigWig file in Galaxy" />
</a>

Click on "main" and the UCSC browser should open within a new window, displaying the data set that you chose.
The default setting for bigWig files is the "dense" display that looks like a heatmap.

<a href="https://raw.github.com/fidelram/deepTools/master/examples/Gal_FAQ_UCSC01.png" target="_blank">
     <img src="https://raw.github.com/fidelram/deepTools/master/examples/Gal_FAQ_UCSC01.png" Title="Screenshot of UCSC browser display of bigWig files" />
</a>

If you would like to display the continuous profile in a "valley-mountain" fashion like the one shown in the IGV screenshot, go to the drop-down menu underneath your custom track and choose "full".

UCSC has large amounts of public data that you can display which you can find by scrolling down the page, beyond your custom track entry. For more information on how to use the UCSC Genome Browser, go [here](https://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html "UCSC user guide").

**Known issues with UCSC**

* **chromosome naming**: UCSC expects chromosome names to be indicated in the format "chr"Number, e.g. chr1. If you mapped your reads to a non-UCSC-standard genome, chances are that chromosomes are labeled just with their number. bigWig files generated from these BAM files will not be recognized by UCSC, i.e. you will see the data set name, but no signal.
* **no upload of bigWig files from your hard drive**: to minimize the computational strains, UCSC relies on streaming bigWig files (i.e. there's no need to load the entire file at once, the browser will always just load the data for the specific region a user is looking at).

What's the best way to integrate the deepTools results with other downstream analyses (outside of Galaxy)?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* you can save all the data tables underlying every image produced by deepTools, i.e. if you would like to plot the average profiles in a different way, you could download the corresponding data (after ticking the profiler option at "advanced output options") and import them into R, Excel, GraphPadPrism etc.


How can I determine basic parameters of a BAM file, such as the number of reads, [read][] length, duplication rate and average DNA fragment length?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Even though [MACS][] is meant to do peak calling for you, it also outputs a number of useful information such as those listed above.
Simply run MACS on the BAM file that you would like to gain the information for and check the .xls file from the MACS output. It will list:

* tag length = read length
* duplication rate
* number of tags = number of reads
* d = distance = average DNA fragment size

--------------------------------------------
[bamCorrelate]: https://github.com/fidelram/deepTools/wiki/QC
[bamFingerprint]: https://github.com/fidelram/deepTools/wiki/QC
[computeGCBias]: https://github.com/fidelram/deepTools/wiki/QC
[bamCoverage]: https://github.com/fidelram/deepTools/wiki/Normalizations
[bamCompare]: https://github.com/fidelram/deepTools/wiki/Normalizations
[correctGCbias]: https://github.com/fidelram/deepTools/wiki/Normalizations
[computeMatrix]: https://github.com/fidelram/deepTools/wiki/Visualizations
[heatmapper]: https://github.com/fidelram/deepTools/wiki/Visualizations
[profiler]: https://github.com/fidelram/deepTools/wiki/Visualizations
[MACS]: http://www.ncbi.nlm.nih.gov/pubmed/22936215 "How to use MACS, Nature Protocols"
[CCAT]: http://www.ncbi.nlm.nih.gov/pubmed/20371496 "CCAT original publication"
[SICER]: http://bioinformatics.oxfordjournals.org/content/25/15/1952.full "SICER original publication"

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

