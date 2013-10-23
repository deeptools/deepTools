======================================================================
deepTools
======================================================================
### user-friendly tools for the normalization and visualization of deep-sequencing data


deepTools addresses the challenge of handling the large amounts of data 
that are now routinely generated from sequencing centers. To do so, deepTools contains useful modules to process the mapped
reads data to create coverage files in standard bedGraph and bigWig file formats. By doing so, deepTools allows the creation of normalized coverage files or the comparison between two files (for example, treatment and control). Finally, using such normalized and standardized files, multiple
visualizations can be created to identify enrichments with
functional annotations of the genome. For a gallery of images that
can be produced, see
http://f1000.com/posters/browse/summary/1094053

For support, questions, or feature requests contact: deeptools@googlegroups.com

### Table of Contents  
[How we use deepTools](#weUse)
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

[Glossary](https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing)

<a name="weUse"/>
How we use deepTools
--------------------------------
The majority of samples that we handle within our facility come from ChIP-seq experiments, therefore you will find many examples from ChIP-seq analyses. This does not mean that deepTools is restricted to ChIP-seq data analysis, but some tools, such as _bamFingerprint_ specifically address ChIP-seq-issues. (That being said, we do process quite a bit of RNA-seq and genomic sequencing data, too.)

Our work usually begins with one or more [FASTQ][] file(s) of deeply-sequenced samples. After a first quality control using [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ "Check out FASTQC"), we align the reads to the reference genome, e.g. using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml "bowtie, one of the most popular aligners").
We then use deepTools to assess the quality of the aligned reads:

1. __Correlation between [BAM][] files__ (_bamCorrelate_). This is a very basic test to see whether the sequenced and aligned reads meet your expectations. We use this check to assess the reproducibility - either between replicates and/or between different experiments that might have used the same antibody/the same cell type etc. For instance, replicates should correlate better than differently treated samples.
2. __GC bias check__ (_computeGCbias_). 
3. __ChIP strength__

Once we're satisfied by the basic quality checks, we normally convert the large [BAM][] files into a leaner data format, typically [bigWig][]. bigWig files have many advantages over BAM files:
+ smaller in size
  - useful for data sharing & storage
  - intuitive visualization in Genome Browsers (e.g. UCSC Genome Browser, IGV)
  - more efficient downstream analyses are possible

The deepTools modules _bamCompare_ and _bamCoverage_ do not only allow the simple conversion from BAM to bigWig (or [bedGraph][] for that matter), the main reason why we developed those tools was that we wanted to be able to __normalize__ the read coverages so that we could compare different samples despite differences in sequencing depth, GC biases and so on.

Finally, once all the files have passed our visual inspections, the fun aka downstream analyses with _heatmapper_ and _profiler_ can begin! 

 
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

 

Here's a concise summary of the tools:

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


------------------------------------
[BAM]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "binary version of a SAM file; contains all information about aligned reads"
[SAM]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file containing all information about aligned reads"
[bigWig]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "binary version of a bedGraph file; contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[bedGraph]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file that contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[FASTQ]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file of raw reads (almost straight out of the sequencer)"
### References
[Benjamini and Speed]: http://nar.oxfordjournals.org/content/40/10/e72 "Nucleic Acids Research (2012)"
[Diaz et al.]: http://www.degruyter.com/view/j/sagmb.2012.11.issue-3/1544-6115.1750/1544-6115.1750.xml "Stat. Appl. Gen. Mol. Biol. (2012)"


This tool is developed by the [Bioinformatics Facility](http://www1.ie-freiburg.mpg.de/bioinformaticsfac) at the [Max Planck Institute for Immunobiology and Epigenetics, Freiburg](http://www1.ie-freiburg.mpg.de/).
