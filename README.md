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

![gallery](https://raw.github.com/fidelram/deepTools/master/examples/collage.png)

### Table of Contents  

[What can I do with deepTools? Overview!](#usage)

[How to install deepTools](#installation) 

[How we use deepTools](#weUse)

[Basic options and parameters of deepTools](#parameters)  



More detailed information about the individual programs:
  + [deepTools for data quality checks](https://github.com/fidelram/deepTools/blob/master/manual/QC.md)
  + [deepTools for normalizations](https://github.com/fidelram/deepTools/blob/master/manual/normalizations.md)
  + [deepTools for visualizations](https://github.com/fidelram/deepTools/blob/master/manual/visualizations.md)

[FAQ](#FAQ)

[Glossary](https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing)

---------------------------------------------------------------------------------------------------------------------

<a name="usage"/>
What can I do with deepTools?
---------------

deepTools consists of a set of modules that can be used independently to work with mapped reads. We have subdivided such tasks into *quality controls*, *normalizations* and *visualizations*.

Here's a concise summary of the tools - if you would like more detailed information about the individual tools and example figures, follow the links in the table or check out our general description of [how we use deepTools](#weUse)

| tool | type | input files | main output file(s) | application |
|------|------|-------------|---------------------|-------------|
| [bamCorrelate][] | QC | 2 or more BAM | clustered heatmap | Pearson or Spearman correlation between read distributions |
| [bamFingerprint][] | QC | 2 BAM | 1 diagnostic plot | assess enrichment strength of a ChIP sample |
| [computeGCBias][] | QC | 1 BAM | 2 diagnostic plots | calculate the exp. and obs. GC distribution of reads|
| [bamCoverage][] | normalization | BAM | bedGraph or bigWig | obtain the normalized read coverage of a single BAM file |
| [bamCompare][] | normalization | 2 BAM | bedGraph or bigWig | normalize 2 BAM files to each other using a mathematical operation of your choice (e.g. log2ratio, difference)|
| [computeMatrix][] | visualization | 1 bigWig, 1 BED | zipped file, to be used with heatmapper or profiler | compute the values needed for heatmaps and summary plots |
| [heatmapper][] | visualization | computeMatrix output | heatmap of read coverages | visualize the read coverages for genomic regions |
| [profiler][] | visualization | computeMatrix output | summary plot ("meta-profile") | visualize the average read coverages over a group of genomic regions |


<a name="installation"/>
Installation
---------------
[Installation from source](#linux)

[Installation on a Mac](#mac)

[Troubleshooting](#trouble)

[Installation via Galaxy](#galaxy)


<a name="linux"/>
### Installation from source (Linux, command line)

The easiest way to install deepTools is by __downloading the source file and using python pip__ or easy_install tools:

Requirements: Python 2.7, numpy, scipy installed

Commands:

      $ cd ~
      $ export PYTHONPATH=$PYTHONPATH:~/lib/python2.7/site-packages
      $ export PATH=$PATH:~/bin:~/.local/bin

If pip is not already available, install with:

      $ easy_install --prefix=~ pip

Install deepTools and dependencies with pip:

      $ pip install --user deeptools
Done.




__Another option is to clone the repository:__
	
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

<a name="mac"/>
### Installation on a MAC

Requirement: Python 2.7 installed

Numpy, Scipy, matplotlib are also required - in case you haven't installed them yet, download the packages and install them using dmg images:
- http://sourceforge.net/projects/numpy/files/NumPy/
- http://sourceforge.net/projects/scipy/files/scipy/
- http://matplotlib.org/downloads.html

Then install deepTools via the terminal ("Applications" --> "Terminal"):

     $ cd ~
     $ export PYTHONPATH=$PYTHONPATH:~/lib/python2.7/site-packages
     $ export PATH=$PATH:~/bin:~/.local/bin:~/Library/Python/2.7/bin

If pip is not already available, install with:

     $ easy_install --prefix=~ pip

Install deepTools and dependencies with pip:

     $ pip install --user pysam
     $ pip install --user bx-python
     $ pip install --user --no-deps deeptools


<a name="trouble"/>
##### Troubleshooting
The easy_install command is provided by the python package setuptools.
You can download the package from https://pypi.python.org/pypi/setuptools

     $ wget https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py -O - | python
     
or the user-specific way:

     $ wget https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py
     $ python ez_setup.py --user

Numpy/Scipy Installation:
http://www.scipy.org/install.html

<a name="galaxy"/>
#### Galaxy Installation

deepTools can be easily integrated into [Galaxy](http://galaxyproject.org). All wrappers and dependencies are 
available in the [Galaxy Tool Shed](http://toolshed.g2.bx.psu.edu/view/bgruening/deeptools).


##### Installation via Galaxy API (recommended)

At first generate an [API Key](http://wiki.galaxyproject.org/Admin/API#Generate_the_Admin_Account_API_Key) for your admin 
user and run the the installation script:

	python ./scripts/api/install_tool_shed_repositories.py --api YOUR_API_KEY -l http://localhost:8080 --url http://toolshed.g2.bx.psu.edu/ -o bgruening -r <revision> --name deeptools --tool-deps --repository-deps --panel-section-name deepTools

The -r argument specifies the version of deepTools. You can get the latest revsion number from the test tool shed or with the following command:

	hg identify http://toolshed.g2.bx.psu.edu/view/bgruening/deeptools

You can watch the installation status under: Top Panel → Admin → Manage installed tool shed repositories


##### Installation via webbrowser

- go to the [admin page](http://localhost:8080/admin)
- select *Search and browse tool sheds*
- Galaxy tool shed → Sequence Analysis → deeptools
- install deeptools

remember: for support, questions, or feature requests contact: deeptools@googlegroups.com

<a name="weUse"/>
How we use deepTools
--------------------------------
The majority of samples that we handle within our facility come from ChIP-seq experiments, therefore you will find many examples from ChIP-seq analyses. This does not mean that deepTools is restricted to ChIP-seq data analysis, but some tools, such as _bamFingerprint_ specifically address ChIP-seq-issues. (That being said, we do process quite a bit of RNA-seq, other -seq and genomic sequencing data using deepTools, too.)

[Here](https://docs.google.com/file/d/0B8DPnFM4SLr2UjdYNkQ0dElEMm8/edit?usp=sharing "From aligned reads to coverage profiles using deepTools") are slides that we used for teaching at the University of Freiburg.

As depicted in the figure down below, our work usually begins with one or more [FASTQ][] file(s) of deeply-sequenced samples. After a first quality control using [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ "Check out FASTQC"), we align the reads to the reference genome, e.g. using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml "bowtie, one of the most popular aligners").
We then use deepTools to assess the quality of the aligned reads:

1. __Correlation between [BAM][] files__ (_bamCorrelate_). This is a very basic test to see whether the sequenced and aligned reads meet your expectations. We use this check to assess the reproducibility - either between replicates and/or between different experiments that might have used the same antibody/the same cell type etc. For instance, replicates should correlate better than differently treated samples.
2. __GC bias check__ (_computeGCbias_). Many sequencing protocols require several rounds of PCR-based amplification of the DNA to be sequenced. Unfortunately, most DNA polymerases used for PCR introduce significant GC biases as they prefer to amplify GC-rich templates. Depending on the sample (preparation), the GC bias can vary significantly and we routinely check its extent. In case we need to compare files with different GC biases, we use the _correctGCbias_ module to match the GC bias.
See the paper by [Benjamini and Speed][] for many insights into this problem.
3. __Assessing the ChIP strength__. This is a QC we do to get a feeling for the signal-to-noise ratio in samples from ChIP-seq experiments. It is based on the insights published by [Diaz et al.][].

Once we're satisfied by the basic quality checks, we normally __convert the large [BAM][] files into a leaner data format, typically [bigWig][]__. bigWig files have several advantages over BAM files that mainly stem from their significantly decreased size:
  - useful for data sharing & storage
  - intuitive visualization in Genome Browsers (e.g. UCSC Genome Browser, IGV)
  - more efficient downstream analyses are possible

The deepTools modules _bamCompare_ and _bamCoverage_ do not only allow the simple conversion from BAM to bigWig (or [bedGraph][] for that matter), __the main reason why we developed those tools was that we wanted to be able to *normalize* the read coverages__ so that we could compare different samples despite differences in sequencing depth, GC biases and so on.

Finally, once all the files have passed our visual inspections, the fun of downstream analyses with _heatmapper_ and _profiler_ can begin! 

Here's a visual summary of our average workflow - deepTools modules are indicated in bold letters, alternative software such as FASTQC and bowtie are noted in regular font. Everything written in red is related to quality control (QC) of the samples.

![flowChartI](https://raw.github.com/fidelram/deepTools/master/examples/flowChart_BAMtoBIGWIG.png "Average analysis and QC workflow")

 

<a name="parameters"/>
General information about deepTools usage
---------------------------------------------------------

All tools require the user to specify input files, output file names, optional and mandatory parameters.

Here we point out some parameters that you might find especially useful in your regular usage of deepTools.

#### Parameters to decrease the run time
  + numberOfProcessors 
  + region - in case you're testing whether a certain plot works and gives you the output you're hoping for, you can speed things up by focusing on a certain genome region, e.g. chr4 or chr2:100000200000

#### filtering BAMs while processing
  + ignoreDuplicates
  + minMappingQuality

#### random tips 
  + output format of plots should be indicted by the file ending, e.g. MyPlot.pdf will return a pdf, MyPlot.png a png-file
  + all tools that produce plots can also output the underlying data - this can be useful in case you don’t like the deepTools visualization as you can then use the data matrices produced by deepTools with your favourite plotting module, e.g. R or Excel


<a name="FAQ"/>
FAQs
-------
#### How does deepTools handle data from paired-end sequencing?
Generally, all the modules working with BAM files (_bamCorrelate, bamCoverage, bamCompare, bamFingerprint, computeGCbias_)
recognize paired-end sequencing data. You can enforce to ignore the fragment length based on the mate pairs using the option __doNotExtendPairedEnds_

#### Where can I download the 2bit genome files required for _computeGCbias_?
The 2bit files of most genomes can be found [here](http://hgdownload.cse.ucsc.edu/gbdb/).
Search for the .2bit ending. Otherwise, __fasta files can be converted to 2bit__ using the UCSC programm
faToTwoBit (available for different plattforms from [here](http://hgdownload.cse.ucsc.edu/admin/exe/)

#### When should I exclude peaks from the GC bias computation?
Here's what we do: We usually check the GC bias (using _computeGCbias_) on the entire data set. If we notice that the majority of the genome is dramatically biased (i.e. instead of a straight line for genome regions of 30-60% GC content we see a diagonale or something even weirder) and we want to compare this particular, biased sample with another unbiased (or differently GC biased) sample, then we need to correct the GC bias (using _correctGCbias_).

Now, imagine that the biased sample is a ChIP for a protein binding to methylated or mammalian promoter regions. By default, such data is going to have a GC bias because CpG-rich regions should be enriched. Therefore, you can use the option to exclude regions of signficant enrichment from the GC bias computation and subsequent correction to only account for the background bias without interfering with the true signal too much. 

------------------------------------
[BAM]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "binary version of a SAM file; contains all information about aligned reads"
[SAM]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file containing all information about aligned reads"
[bigWig]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "binary version of a bedGraph file; contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[bedGraph]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file that contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[FASTQ]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file of raw reads (almost straight out of the sequencer)"

[bamCorrelate]: https://github.com/fidelram/deepTools/blob/master/manual/QC.md
[bamFingerprint]: https://github.com/fidelram/deepTools/blob/master/manual/QC.md
[computeGCBias]: https://github.com/fidelram/deepTools/blob/master/manual/QC.md

[bamCoverage]: https://github.com/fidelram/deepTools/blob/master/manual/normalizations.md
[bamCompare]: https://github.com/fidelram/deepTools/blob/master/manual/normalizations.md

[computeMatrix]: https://github.com/fidelram/deepTools/blob/master/manual/visualizations.md
[heatmapper]: https://github.com/fidelram/deepTools/blob/master/manual/visualizations.md
[profiler]: https://github.com/fidelram/deepTools/blob/master/manual/visualizations.md
### References
[Benjamini and Speed]: http://nar.oxfordjournals.org/content/40/10/e72 "Nucleic Acids Research (2012)"
[Diaz et al.]: http://www.degruyter.com/view/j/sagmb.2012.11.issue-3/1544-6115.1750/1544-6115.1750.xml "Stat. Appl. Gen. Mol. Biol. (2012)"


This tool suite is developed by the [Bioinformatics Facility](http://www1.ie-freiburg.mpg.de/bioinformaticsfac) at the [Max Planck Institute for Immunobiology and Epigenetics, Freiburg](http://www1.ie-freiburg.mpg.de/).
