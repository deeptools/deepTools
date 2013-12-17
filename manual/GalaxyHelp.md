Using deepTools within Galaxy
=============================

We have a publicly available deepTools installation embedded within the Galaxy framework: __deeptools.ie-freiburg.mpg.de__

[Galaxy][] is a tremendously useful platform developed by the Galaxy Team at Penn State.
This platform is meant to offer access to a large variety of bioinformatics tools that __can be used without computer programming experiences__.

We have compiled several information that will hopefully get you started with your data analysis quickly. If you have never worked with Galaxy before, we recommend to start from the top.

If you just need a refresher of __how to upload data into Galaxy__, please have a look [down below](#upload).

If you do not know the difference between a BAM and a BED file, that's fine - just make sure you have a look at this brief overview [here](https://drive.google.com/file/d/0B8DPnFM4SLr2UHY3cHNZdTFEcDg/edit?usp=sharing "NGS Data Formats") before starting your analysis as high-throughput sequencing data relies on several specific __data formats__.

If you would like to dive right into the analysis of BAM or bigWig files using [deepTools Galaxy][], [here](https://drive.google.com/file/d/0B8DPnFM4SLr2M0lZYkl6UEw2WGs/edit?usp=sharing "Visual Manual of deepTools Galaxy") is a __manual that covers the basic functions of deepTools__, starting from an overview of a typical workflow of NGS data analysis and ending with several different examples to demonstrate the power of heatmaps and summary plots.

-----------------------------------------------------------------------------------
## Table of Content

  * [Basic features of Galaxy](#basics)
  * [Data import](#upload)
      * [upload files](#dataup)
      * [import shared files](#dataim)
      * [download annotation and publicly available tracks](#downloadann)
  * [Tools](#tools)
  * [Example Workflows](#workflows)
  * [deepTools Galaxy Tipps and FAQ](https://github.com/fidelram/deepTools/blob/master/manual/GalaxyFAQs.md)
  * [Where to get help](#help)
   

<a name="basics">
Basic features of Galaxy
-------------------------

The Galaxy team develops the platform, but since it is impossible to meet all bioinformatics needs (that can range from evolutionary analysis to mass spec data to high-throughput DNA sequencing (and way beyond)) with one single web server, many institutes have installed their own versions of Galaxy tuned to their specific needs. Our deepTools Galaxy is such a specialized server dedicated to the analysis of high-throughput DNA sequencing data. The overall makeup of this web server, however, is the same as for any other Galaxy installation, so if you've used Galaxy before, this section will probably not give you any new insights. 

### The start site

Here is a screenshot of what you should see at deeptools.ie-freiburg.mpg.de:

![startSite](https://raw.github.com/fidelram/deepTools/master/examples/Gal_startsite.png "Screenshot of deepTools Galaxy")

The start site contains 4 main features:

  * __Top menu__: will lead you to other sections of Galaxy (away from the actual analysis part), such as workflows (registered users only) and content shared with you by other users such as sample data sets, pages and workflows
  * __Tool panel__ "What can be done": via this menu you can find all the _tools_ installed in this Galaxy instance
  * __Main frame__ "What am I doing?": the center frame is your main working space where input will be required from you once you use a tool. In addition, you will always find general information about the tool here
  * __History panel__ "What did I do?": here you can find all _files_ that one produces or uploads
      * the history is like a log book: everything you ever did is recorded here (unless you deleted things permanently)
      * histories can be shared with other users, they can also be downloaded
      * for each file that was produced, you will find all kinds of useful information such as the tool that was used to create the file, the tool's parameters etc.

For those visual learners, here's an annotated screenshot:
        
![startSiteAnnotated](https://raw.github.com/fidelram/deepTools/master/examples/Gal_startsite_with_comments.png "Screenshot of deepTools Galaxy")

In the default state of the tool panel you see the __tool categories__, e.g. "Get Data". If you click on them, you will see the __individual tools__ belonging to each category, e.g. "Upload File from your computer", "UCSC Main table browser" and "Biomart central server" in case you clicked on "Get Data". To use a tool such as "Upload File from your computer", just click on it.

The tool search panel is extremely useful as it allows you to enter a  key word (e.g. "BAM") that will lead to all the tools mentioning the key word in their tool name.

Once you've uploaded any kind of data, you will find the history on the right hand side filling up with green tiles. Each tile corresponds to one data set that you either uploaded or created. The data sets can be images, raw sequencing files, text files, tables, virtually anything. The content of each data set cannot be modified - everytime you want to change something _within_ a data file (e.g. you would like to sort the values or add a line or cut a column), you will have to use a Galaxy tool that will lead to a _new_ data set being produced. Every data set can be downloaded to your computer. 

Have a look at the following screenshot to get a feeling for how many information Galaxy keeps for you (which makes it very feasible to reproduce any given analysis):

![DataSet](https://raw.github.com/fidelram/deepTools/master/examples/Gal_screenshot_dataSet.png "Galaxy data set")


Each data set can have 4 different states that are intuitively color-coded:
![DataSetStates](https://raw.github.com/fidelram/deepTools/master/examples/Gal_screenshot_dataSetStates.png "States of Galaxy data sets")


If you encounter a failure after you've run a tool, please follow those steps (in this order):

1. click on the center button on the lower left corner of the failed data set ("i"): now check whether you chose the __correct data files__
2. if you're sure that you chose the correct files, hit the re-run button (blue arrow in the lower left corner) - check again whether your files had the __correct file format__ 
  + if you suspect that the format might be incorrectly assigned (e.g. a file that should be a bed-file is labelled as a tabular file), click the edit button of the input data file - ther you can change the corresponding attributes
3. if you've checked your input data and the error is persisting, click on the green bug (lower left corner of the failed data set) and send the __bug report__ to us.


<a name="upload">
Data import into Galaxy
-------------------------

There are three main ways to populate your Galaxy history with data files:

1. [Data upload from your computer](#dataup)
2. [Import a shared data set from the Galaxy data library](#dataim)
3. [Download annotation data from public servers](#downloadann)
4. [__For registered users only__: Copy data sets between histories](#copy)


<a name="dataup">
#### Upload files from your computer
The data upload of files <2 GB that lie on your computer is fairly straight-forward: click on the category "Get data" and choose the tool "Upload file". Then select the file via the "Browse" button.
        
![DataUpload](https://raw.github.com/fidelram/deepTools/master/examples/Gal_DataUpload.png "Screenshot of the data upload tool")


For files >2GB there's the option to upload via an FTP server. If your data is available via an URL that links to an FTP server, you can simply paste the URL in the empty text box. 

If you do not have access to an FTP server, you can directly upload to our Galaxy's FTP. For this, you must first register with deeptools.ie-freiburg.mpg.de (via “User” --> “register”; registration requires an email address and is free of charge). You will also need an FTP client, e.g. [filezilla](https://filezilla-project.org/ "Filezilla web site"). Once you've registered with us and you've got the FTP client, login to the __FTP client__ using your deepTools Galaxy user name and password (host: deeptools.ie-freiburg.mpg.de). Down below you see a screenshot of how that looks like with filezilla. Move the file you wish to upload to the remote site and go back to [deepTools Galaxy][] where the file should now appear once you click on the tool "Upload file" (--> "Files uploaded via FTP"). Hit “execute”.

![Filezilla](https://raw.github.com/fidelram/deepTools/master/examples/Gal_filezilla.png "Screenshot of filezilla")


<a name="dataim">
#### Import data sets from the Galaxy data library

If you would like to try out some tools and play around with sample data, you can import some files that we have saved within the general data storage of the deepTools Galaxy server. These files are not part of anyone's history in particular, but everyone can import them into his or her own history.

You can reach the data library via "Shared Data" in the top menu where you select "Data Libraries". (As you can see, there's more than just data that can be shared within Galaxy.)

Within the Data Library you will find a folder called "Sample Data" that contains data that we downloaded from the [Roadmap project][]. More precisely, we donwloaded the [FASTQ][] files and mapped the reads to the human reference genome (version hg19). The [BAM] files

![DataLibrary](https://raw.github.com/fidelram/deepTools/master/examples/Gal_DataLib.png "Screenshots of how to get to the Data Library")


<a name="downloadann">
#### Download annotation files from public data bases

In many cases you will want to query your sequencing data results for known genome annotation, such as genes, exons, transcription start sites etc. These information can be obtained via the two main sources of genome annotation, [UCSC][] and [BioMart][]. Please note that UCSC and BioMart will cater to different ways of genome annotation, i.e. genes defined in UCSC might not correspond to the same regions in a gene file downloaded from BioMart. (For a brief overview over the issues of genome annotation, you can check out [Wikipedia](http://en.wikipedia.org/wiki/Genome_project "Genome annotation article"), if you'd always wanted to know much more about those issues, [this](http://www.ncbi.nlm.nih.gov/pubmed/22510764 "A beginner's guide to eukaryotic genome annotation, Nat. Genetics, 2012") might be a good start.)

You can access the data stored at UCSC or BioMart conveniently through our Galaxy instance which will import the resulting files into your history. Just go to __"Get data"__ --> "UCSC" or "BioMart".

The majority of annotation files will probably be in [BED][] format, however, there is no limitation on the kind of data you would like to incorporate. UCSC, for example, offers a wide range of data that you can browse via the "group" and "track" menus (for example, you could download the GC content of the genome as a signal file from UCSC via the "group" menu ("Mapping and Sequencing Tracks")). 

So, here's a screenshot from downloading a BED-file of all RefSeq genes defined for the human genome (version hg19):
![UCSC](https://raw.github.com/fidelram/deepTools/master/examples/Gal_UCSC.png "Screenshot of the UCSC main table browser")


And here's how you would do it for the BioMart approach:
![Biomart](https://raw.github.com/fidelram/deepTools/master/examples/Gal_biomart.png "Screenshot of the Biomart main table browser")

Per default, BioMart will not output a BED file like UCSC does. It is therefore important that you make sure you get all the information you need (most likely: chromosome, gene start, gene end, ID, strand) via the "Attributes" section. You can click on the "Results" button at any time to check the format of the table that will be sent to Galaxy (Note that the strand information will be decoded as 1 for "forward" or "plus" strand and -1 for "reverse" or "minus" strand.)

>be aware, that BED files from UCSC will have chromosomes labelled with “chr” while ENSEMBL usually returns just the number – this might lead to incompatibilities, i.e. when working with annotations from UCSC and ENSEMBL, you need to make sure to use the same naming!


<a name="copy">
#### __For registered users only__: Copy data sets between histories
In case you have registered with deepTools Galaxy you can have more than one history. In order to minimize the disk space you're occupying we strongly suggest to __copy__ data sets between histories in case that you would like to use one data set in different histories. This can easily be done via the History panel's option button --> "Copy dataset". In the main frame, you should now be able to select the history you would like to copy from on the left hand side while you should indicate the target history on the right hand side.


<a name="tools">
Which tools can I find in the deepTools Galaxy?
-------------------------------------------------------

As mentioned above, each Galaxy installation can be tuned to the individual interests. Our goal is to provide a Galaxy that enables you to __quality check, process and normalize and subsequently visualize your data obtained by high-throughput DNA sequencing__.

* [deepTools - NGS data handling](#deepTools)
* [peak calling (ChIP-seq specific)](#peaks)
* [operating on genomic intervals](#BED)
* [working with text files and tables](#textfiles)


<a name="deepTools">
#### deepTools

The most important category is called __"deepTools"__ that contains 8 major tools:

| tool | type | input files | main output file(s) | application |
|------|------|-------------|---------------------|-------------|
| bamCorrelate | QC | 2 or more BAM | clustered heatmap | Pearson or Spearman correlation between read distributions |
| bamFingerprint | QC | 2 BAM | 1 diagnostic plot | assess enrichment strength of a ChIP sample |
| computeGCBias | QC | 1 BAM | 2 diagnostic plots | calculate the exp. and obs. GC distribution of reads|
| bamCoverage | normalization | BAM | bedGraph or bigWig | obtain the normalized read coverage of a single BAM file |
| bamCompare | normalization | 2 BAM | bedGraph or bigWig | normalize 2 BAM files to each other using a mathematical operation of your choice (e.g. log2ratio, difference)|
| computeMatrix | visualization | 1 bigWig, 1 BED | zipped file, to be used with heatmapper or profiler | compute the values needed for heatmaps and summary plots |
| heatmapper | visualization | computeMatrix output | heatmap of read coverages | visualize the read coverages for genomic regions |
| profiler | visualization | computeMatrix output | summary plot ("meta-profile") | visualize the average read coverages over a group of genomic regions |


If you would like to have more details about the individual tools, [here](https://drive.google.com/file/d/0B8DPnFM4SLr2M0lZYkl6UEw2WGs/edit?usp=sharing "Visual Manual of deepTools Galaxy") is a _manual that covers the basic functions of deepTools_, starting from an overview of a typical workflow of NGS data analysis and ending with several different examples to demonstrate the power of heatmaps and summary plots.

For each tool, there will specific explanations within the [deepTools Galaxy][] main frame, too.


<a name="peaks">
#### Peak calling

In ChIP-seq analysis, peak calling algorithms are essential downstream analysis tools to identify regions of significant enrichments (i.e. where the ChIP sample contained significantly more sequenced reads than the input control sample). By now, there must be close to 100 programs out there (see [Wilbanks et al.](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0011471 "Wilbanks et al.") for a comparison of peak calling programs).

In contrast to deepTools that were developed for handling and generating _continuous_ genome-wide profiles, peak calling will result in a _list of genomic regions_. Have a look at the screenshot to understand the difference.

![PeaksBigWigs](https://raw.github.com/fidelram/deepTools/master/examples/Gal_peaksVsBigWigs.png "Screenshot of the Biomart main table browser")

We have included the peak callers [MACS][], [SICER][] and [CCAT][] within our Galaxy instance with [MACS][] being the most popular peak calling algorithm for the identification of localized transcription factor binding sites while [SICER][] was developed for diffuse ChIP-seq signals.


<a name="BED">
#### Working with genomic intervals

Galaxy has 2 file formats to store lists of genomic regions:

* INTERVAL
    * tab-separated
    * requirements:
        1. Column: chromosome
        2. Column: start position
        3. Column: end position
     * all other columns can contain any value or character

* BED
    * very similar to INTERVAL, but stricter when it comes to what is expected to be kept in which column:
        * 1. to 3. Column: same as interval
        * Column 4: name
        * Column 5: score
        * Column 6: strand

In case you would like to work with several lists of genomic regions, e.g. generate a new list of regions that are found in two different files etc., there are two categories of tools dedicated to performing these tasks:
* Operate on genomic intervals
* BEDtools

Each tool's function is explained within Galaxy. Do browse those tools as they will give you a very good glimpse of the scope of possible analyses!


<a name="textfiles">
#### Working with text files and tables
In addition to deepTools that were specifically developed for the handling of NGS data, we have incorporated several standard Galaxy tools that enable you to manipulate tab-separated files such as gene lists, peak lists, data matrices etc.

There are 3 main categories:
* __Text manipulation__
    * unlike Excel where you can easily interact with your text and tables via the mouse, data manipulations within Galaxy are strictly based on commands. If you feel like you would like to do something to certain _columns_ of a data set, go through the tools of this category
    * e.g. adding columns, cutting columns, pasting two files side by side, selecting random lines etc. 
* __Filter and Sort__
    * in addition to the common sorting and filtering, there's the very useful toool to __select lines that match an expression"
* __Join, Subtract, Group__
   * this category is very useful if you have several data sets that you would like to work with, e.g. by comparing them


<a name="help">
Example workflows
--------------------
Workflows are Galaxy's equivalent of protocols. This is an extremely useful feature as it allows users to share their protocols and bioinformatic analyses in a very easy and transparent way.
This is the graphical representation of a Galaxy workflow that can easily be modified via drag'n'drop within the workflows manual (you must be registered with deepTools Galaxy to be able to generate your own workflows).
![workflow](https://raw.github.com/fidelram/deepTools/master/examples/Gal_workflow.png "Exemplary Galaxy workflow")

We have compiled several workflows that should give you a feeling for the kinds of analyses that can be done using [deepTools Galaxy].


<a name="help">
Where to get help?
--------------------

Please check our [deepTools Galaxy FAQs](https://github.com/fidelram/deepTools/blob/master/manual/GalaxyFAQs.md)

* general Galaxy help: wiki.galaxyproject.org/Learn
* specific help with deepTools Galaxy: deeptools@googlegroups.com
* if you encounter a failing data set, please send a bug report via Galaxy and we will get in touch


-----------------------------------------------------------------------------------

This tool is developed by the [Bioinformatics Facility](http://www1.ie-freiburg.mpg.de/bioinformaticsfac) at the [Max Planck Institute for Immunobiology and Epigenetics, Freiburg](http://www1.ie-freiburg.mpg.de/).


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
