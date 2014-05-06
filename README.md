======================================================================
deepTools
======================================================================
### user-friendly tools for the normalization and visualization of deep-sequencing data


deepTools addresses the challenge of handling the large amounts of data 
that are now routinely generated from DNA sequencing centers. To do so, deepTools contains useful modules to process the mapped reads data to create coverage files in standard bedGraph and bigWig file formats. By doing so, deepTools allows the creation of **normalized coverage files** or the comparison between two files (for example, treatment and control). Finally, using such normalized and standardized files, multiple
**visualizations** can be created to identify enrichments with
functional annotations of the genome.

For support, questions, or feature requests contact: deeptools@googlegroups.com

Citation: Fidel Ramírez, Friederike Dündar, Sarah Diehl, Björn A. Grüning, and Thomas Manke. [_deepTools: a flexible platform for exploring deep-sequencing data._](http://nar.oxfordjournals.org/content/early/2014/05/05/nar.gku365.abstract) Nucl. Acids Res. first published online May 5, 2014 doi:10.1093/nar/gku365

![gallery](https://raw.github.com/fidelram/deepTools/master/examples/collage.png)

Our [wiki page](https://github.com/fidelram/deepTools/wiki) contains more information on **why we built deepTools**, details on the **individual tool scopes and usages** and an introduction to our deepTools Galaxy web server. It also contains an [FAQ][] section that we update regularly. For more specific troubleshooting, feedback, and tool suggestions, contact us via deeptools@googlegroups.com


-------------------------------------------------------------------------------------------------------------------

<a name="installation"/></a>
Installation
---------------

deepTools are available for:

* command line usage
* integration into Galaxy servers

Details on the installation routines can be found here.

[General Installation](#general)

[Installation on a Mac](#mac)

[Galaxy installation](#galaxy)


<a name="general"/></a>
### General Installation

The easiest way to install deepTools is by using python `pip` or `easy_install tools`:

Requirements: Python 2.7, numpy, scipy (http://www.scipy.org/install.html) installed

Commands:

      $ pip install deeptools
Done.


__A second option is to clone the repository:__
	
	$ git clone https://github.com/fidelram/deepTools
	$ cd deepTools
	$ python setup.py install
	
By default, the script will install the python library and executable
codes globally, which means you need to be root or administrator of
the machine to complete the installation. If you need to
provide a nonstandard install prefix, or any other nonstandard
options, you can provide many command line options to the install
script.

	$ python setup.py --help

For example, to install under a specific location use:

	$ python setup.py install --prefix <target directory>

<a name="mac"></a>
### Installation on a MAC

The easiest way to get numpy and scipy dependencies is to install the
[Anaconda Scientific Python Distribution][]. After installation, open
a terminal ("Applications" → "Terminal") and follow the [General Installation](#general)
  	   
If individual installation of the dependencies is preferred, follow 
those steps:

Requirement: Python 2.7 installed

Download the packages and install them using dmg images:
- http://sourceforge.net/projects/numpy/files/NumPy/
- http://sourceforge.net/projects/scipy/files/scipy/

Then open terminal ("Applications" → "Terminal")
and follow the [General Installation](#general)


<a name="trouble"/></a>

<a name="galaxy"/></a>
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

------------------------------------
[BAM]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "binary version of a SAM file; contains all information about aligned reads"
[SAM]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file containing all information about aligned reads"
[bigWig]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "binary version of a bedGraph file; contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[bedGraph]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file that contains genomic intervals and corresponding scores, e.g. average read numbers per 50 bp"
[FASTQ]: https://docs.google.com/document/d/1Iv9QnuRYWCtV_UCi4xoXxEfmSZYQNyYJPNsFHnvv9C0/edit?usp=sharing "text file of raw reads (almost straight out of the sequencer)"

[bamCorrelate]: https://github.com/fidelram/deepTools/wiki/QC#wiki-bamCorrelate
[bamFingerprint]: https://github.com/fidelram/deepTools/wiki/QC#wiki-bamFingerprint
[computeGCBias]: https://github.com/fidelram/deepTools/wiki/QC#wiki-computeGCbias
[bamCoverage]: https://github.com/fidelram/deepTools/wiki/Normalizations#wiki-bamCoverage
[bamCompare]: https://github.com/fidelram/deepTools/wiki/Normalizations#wiki-bamCompare
[computeMatrix]: https://github.com/fidelram/deepTools/wiki/Visualizations
[heatmapper]: https://github.com/fidelram/deepTools/wiki/Visualizations
[profiler]: https://github.com/fidelram/deepTools/wiki/Visualizations
[FAQ]: https://github.com/fidelram/deepTools/wiki/FAQ

[Benjamini and Speed]: http://nar.oxfordjournals.org/content/40/10/e72 "Nucleic Acids Research (2012)"
[Diaz et al.]: http://www.degruyter.com/view/j/sagmb.2012.11.issue-3/1544-6115.1750/1544-6115.1750.xml "Stat. Appl. Gen. Mol. Biol. (2012)"
[Anaconda Scientific Python Distribution]: https://store.continuum.io/cshop/anaconda/

This tool suite is developed by the [Bioinformatics Facility](http://www1.ie-freiburg.mpg.de/bioinformaticsfac) at the [Max Planck Institute for Immunobiology and Epigenetics, Freiburg](http://www1.ie-freiburg.mpg.de/).

[Wiki Start Page](https://github.com/fidelram/deepTools/wiki) | [deepTools Galaxy](http://deeptools.ie-freiburg.mpg.de) | [FAQ](https://github.com/fidelram/deepTools/wiki/FAQ)
