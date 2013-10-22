======================================================================
deepTools
======================================================================
### user-friendly tools for the normalization and visualization of deep-sequencing data


deepTools addresses the challenge of handling the large amounts of data 
that are now routinely generated from sequencing centers. To do so, deepTools contains useful routines to process the mapped
reads data through removal of duplicates and different filtering options
to create coverage files in standard bedGraph and bigWig file formats.
In addition, deepTools allow the creation of normalized coverage files or the
comparison between two files (for example, treatment and control).
Finally, using such normalized and standardized files, multiple
visualizations can be created to identify enrichments with
functional annotations of the genome. For a gallery of images that
can be produced, see
http://f1000.com/posters/browse/summary/1094053

![example heatmap](https://raw.github.com/fidelram/deepTools/master/examples/heatmaps.png)

For support, questions, or feature requests contact: deeptools@googlegroups.com

### Table of Contents  
[How we typically use deepTools](#weUse)
[How to install deepTools](#installation)  
[Basic options and parameters of deepTools](#parameters)  

[Using deepTools](#usage):

[Data quality checks](#qc)
  * [bamCorrelate](#bamCorrelate)
  * [bamFingerprint](#bamFingerprint)
  * [computeGCbias](#computeGCbias)

[Normalizations](#norm)
  * [correctGCbias](#correctGCbias)
  * [bamCoverage](#bamCoverage)
  * [bamCompare](#bamCompare)

[Visualizations](#visualizations)
  * [computeMatrix](#computeMatrix)
  * [heatmapper](#heatmapper)
  * [profiler](#profiler)

[Glossary](#glossary)

<a name="weUse"/>
How we use deepTools
--------------------------------
Our work usually begins with one or more __FASTQ__ file(s) of deeply-sequenced samples. . After a first quality control using [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), we align the reads to the reference genome, e.g. using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).

 
<a name="installation"/>
Installation
---------------
### Installation from source

The easiest way to install deepTools is by downloading the
source file and using python pip or easy_install tools:

	$ pip install deepTools
	$ vim [deepTools folder]/config/deepTools.cfg

The `deepTools.cfg` file contains several variables that
need to be adjusted.
 
Other option is to clone the repository
	
	$ git clone https://github.com/fidelram/deepTools

Then go to the deepTools directory, edit the `deepTools.cfg` 
file and then run the install script a:

	$ cd deepTools
	$ vim config/deepTools.cfg
	$ python setup.py install
	

By default, the script will install python library and executable
codes globally, which means you need to be root or administrator of
the machine to complete the installation. If you need to
provide a nonstandard install prefix, or any other nonstandard
options, you can provide many command line options to the install
script.

	$ python setup.py --help

To install under a specific location use:

	$ python setup.py install --prefix <target directory>


#### Galaxy Installation with the Tool Shed

deepTools can be easily integrated into [Galaxy](http://galaxyproject.org). All wrappers and dependencies are 
available in the [Galaxy Tool Shed](http://testtoolshed.g2.bx.psu.edu/view/bgruening/deeptools).


#### Installation via Galaxy API (recommended)

At first generate an [API Key](http://wiki.galaxyproject.org/Admin/API#Generate_the_Admin_Account_API_Key) for your admin 
user and run the the installation script:

	python ./scripts/api/install_tool_shed_repositories.py --api YOUR_API_KEY -l http://localhost:8080 --url http://testtoolshed.g2.bx.psu.edu/ -o bgruening -r c8a0dc481493 --name deeptools --tool-deps --repository-deps --panel-section-name deepTools

The -r argument specifies the version of deepTools. You can get the latest revsion number from the test tool shed or with the following command:

	hg identify http://testtoolshed.g2.bx.psu.edu/repos/bgruening/deeptools

You can watch the installation status under: Top Panel → Admin → Manage installed tool shed repositories


#### Installation via webbrowser

- go to the [admin page](http://localhost:8080/admin)
- select *Search and browse tool sheds*
- Galaxy test tool shed → Sequence Analysis → deeptools
- install deeptools

<a name="usage"/>
Using deepTools
---------------

deepTools consists of a set of modules that can be used independently to work with mapped reads. We have subdivided such tasks into *quality controls*, *normalizations* and *visualizations*.

Given, for example, a 

Here's a concise summary of the tools.

| tool | type | input files | main output file(s) | application |
|------|--------|-------------|--------------- |---------------|
| bamCorrelate | QC | 2 or more BAM | clustered heatmap | Pearson or Spearman correlation between read distributions |
| bamFingerprint | QC | 2 BAM | 1 diagnostic plot | assess enrichment strength of a ChIP sample |
| computeGCBias | QC | 1 BAM | 2 diagnostic plots | calculate the exp. and obs. GC distribution of reads|
| bamCoverage | normalisation | BAM | bedGraph or bigWig | obtain the normalized read coverage of a single BAM file |
| bamCompare | normalisation | 2 BAM | bedGraph or bigWig | normalize 2 BAM files to each other using a mathematical operation of your choice (e.g. log2ratio, difference)|
| computeMatrix | visualization | 1 bigWig, 1 BED | zipped file, to be used with heatmapper or profiler | compute the values needed for heatmaps and summary plots |
| heatmapper | visualization | computeMatrix output | heatmap of read coverages | visualize the read coverages for genomic regions |
| profiler | visualziation | computeMatrix output | summary plot ("meta-profile") | visualize the average read coverages over a group of genomic regions |

<a name="parameters"/>
General information about deepTools usage
---------------------------------------------------------

Input files, Output files, optional and mandatory parameters
-	specify the number of processors
-	always specify output names
-	output format of plots should be indicted by the file ending, e.g. MyPlot.pdf will return a pdf, MyPlot.png a png-file
-	all tools that produce plots can also output the underlying data  if you don’t like the visualization deepTools is providing you can always use the data matrices with your favourite plotting module, e.g. R or Excel

<a name="qc"/>
Quality checks
--------------
<a name="bamCorrelate"/>
### bamCorrelate

-	Output files:
o	plot
o	matrix
-	Figure: replicates, input vs ChIP, patient 1 vs patient2

<a name="bamFingerprint"/>
### bamFingerprint
-	What it does (reference to Diaz)
-	Output files: 
o	plot
o	matrix
-	What you can see
-	Figure: good TF vs. bad TF vs. broad histone mark

<a name="computeGCbias"/>
### computeGCbias
-	What it does (reference to Benjamini)
-	Output files
o	plot
o	matrix
-	Figure: strong bias vs. no bias vs. histone mark with inherent GC bias

<a name="norm"/>
Normalizations
---------------------
<a name="correctGCbias"/>
### CorrectGCbias 
-	What it does (uses output from computeGCbias)
-	output files

<a name="bamCoverage"/>
### bamCoverage
-	What it does
-	output files

<a name="bamCompare"/>
### bamCompare
-	What it does
-	output files
-	Figure: IGV snapshots of a) BAM,  b) Norm To 1x, c)difference

<a name="visualizations"/>
Visualizations
--------------

<a name="computeMatrix"/>
### computeMatrix

<a name="heatmapper"/>
### heatmapper

<a name="profiler"/>
### profiler

<a name="glossary"/>
Glossary
------------------------------------
this is to test whether I can see the explanation of a [BAM][] file
and this to test whether I can link several things to the same address: [SAM][]
[BAM]: http://daringfireball.net/projects/markdown/syntax#html "binary version of a SAM file"
[SAM]: http://daringfireball.net/projects/markdown/syntax#html "nonbinary"
### References


This tool is developed by the [Bioinformatics Facility](http://www1.ie-freiburg.mpg.de/bioinformaticsfac) at the [Max Planck Institute for Immunobiology and Epigenetics, Freiburg](http://www1.ie-freiburg.mpg.de/).
