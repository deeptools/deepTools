Using deepTools within Galaxy
=============================

We have a publicly available deepTools installation embedded within the Galaxy framework: __deeptools.ie-freiburg.mpg.de__
[Galaxy][] is a platform developed by the Galaxy Team at Penn State and is an amazing project.
This platform is meant to offer access to a large variety of bioinformatics tools that __can be used without computer programming experiences__.


If you are _new to working with_ Galaxy in general, please have a look at this brief [introduction](https://drive.google.com/file/d/0B8DPnFM4SLr2MzRsb1hPMXBkN00/edit?usp=sharing "Intro into Galaxy").


If you just need a refresher of _how to upload data into Galaxy_, please have a look [here](https://drive.google.com/file/d/0B8DPnFM4SLr2MGI4cHFqVDRTVEE/edit?usp=sharing "Getting Data into Galaxy").

As high-throughput sequencing data relies on several specific _data formats_, we've compiled a short overview [here](https://drive.google.com/file/d/0B8DPnFM4SLr2UHY3cHNZdTFEcDg/edit?usp=sharing "NGS Data Formats").

And finally, [here](https://drive.google.com/file/d/0B8DPnFM4SLr2M0lZYkl6UEw2WGs/edit?usp=sharing "Visual Manual of deepTools Galaxy") is a _manual that covers the basic functions of deepTools_, starting from an overview of a typical workflow of NGS data analysis and ending with several different examples to demonstrate the power of heatmaps and summary plots.

-----------------------------------------------------------------------------------
## Table of Content

  * [Basic features of Galaxy](#basics)
  * [Galaxy tools (very general)](#tools)
  * [Data upload](#upload)
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

Once you've uploaded any kind of data, you will the history on the right hand side filling up with green bricks. Each brick corresponds to one data set that you either uploaded or created. The data sets can be images, raw sequencing files, text files, tables, virtually anything. The content of each data set cannot be modified.Everytime you want to change something _within_ a data file (e.g. you would like to sort the values), you will have to use a Galaxy tool that will lead to a _new_ data set being produced. Every data set can be downloaded do your computer. 

Have a look at the following screenshot to get a feeling for how many information Galaxy keeps for you (which makes it very feasible to reproduce any given analysis):

![DataSet](https://raw.github.com/fidelram/deepTools/master/examples/Gal_screenshot_dataSet.png "Info for an individual Galaxy data set")


Each data set can have 4 different states that are intuitively color-coded:
![DataSetStates](https://raw.github.com/fidelram/deepTools/master/examples/Gal_screenshot_dataSetStates.png "States of Galaxy data sets")


If you encounter a failure after you've run a tool, please follow those steps (in this order):

1. click on the center button on the lower left corner of the failed data set ("i"): now check whether you chose the __correct data files__
2. if you're sure that you chose the correct files, hit the re-run button (blue arrow in the lower left corner) - check again whether your files had the __correct file format__ 
  + if you suspect that the format might be incorrectly assigned (e.g. a file that should be a bed-file is labelled as a tabular file), click the edit button of the input data file - ther you can change the corresponding attributes
3. if you've checked your input data and the error is persisting, click on the green bug (lower left corner of the failed data set) and send the __bug report__ to us.


<a name="tools">
What kinds of tools can I find in the deepTools Galaxy?
-------------------------------------------------------

As mentioned above, each Galaxy installation can be tuned to the individual interests. Our goal is to provide a Galaxy that enables you to __quality check, process and normalize and subsequently visualize your data obtained by high-throughput DNA sequencing__.


<a name="upload">
Data upload into Galaxy
-------------------------

-----------------------------------------------------------------------------------

This tool is developed by the [Bioinformatics Facility](http://www1.ie-freiburg.mpg.de/bioinformaticsfac) at the [Max Planck Institute for Immunobiology and Epigenetics, Freiburg](http://www1.ie-freiburg.mpg.de/).


-----------------------------------------------------------------------------------
[Galaxy]: http://galaxyproject.org/ "General Galaxy platform from Penn State"
[SAM]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file containing all information about aligned reads"
