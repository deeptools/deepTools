======================================================================
deepTools
======================================================================
[![Build Status](https://travis-ci.org/fidelram/deepTools.svg?branch=master)](https://travis-ci.org/fidelram/deepTools) [![Documentation Status](https://readthedocs.org/projects/deeptools/badge/)](http://deeptools.readthedocs.org/) [![PyPI version](https://badge.fury.io/py/deeptools.svg)](https://badge.fury.io/py/deeptools) [![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io)

### user-friendly tools for the normalization and visualization of deep-sequencing data

deepTools addresses the challenge of handling the large amounts of data 
that are now routinely generated from DNA sequencing centers. To do so, deepTools contains useful modules to process the mapped reads data to create coverage files in standard bedGraph and bigWig file formats. By doing so, deepTools allows the creation of **normalized coverage files** or the comparison between two files (for example, treatment and control). Finally, using such normalized and standardized files, multiple
**visualizations** can be created to identify enrichments with
functional annotations of the genome.

For support, questions, or feature requests contact: deeptools@googlegroups.com

Citation: Fidel Ramírez, Friederike Dündar, Sarah Diehl, Björn A. Grüning, and Thomas Manke. [_deepTools: a flexible platform for exploring deep-sequencing data._](http://nar.oxfordjournals.org/content/early/2014/05/05/nar.gku365.abstract) Nucl. Acids Res. first published online May 5, 2014 doi:10.1093/nar/gku365

>Our **[documentation](http://deeptools.readthedocs.org/)** contains more details on the [individual tool scopes and usages](http://deeptools.readthedocs.org/en/latest/content/list_of_tools.html) and an [introduction to our deepTools Galaxy web server](http://deeptools.readthedocs.org/en/latest/content/help_galaxy_intro.html) including [step-by-step protocols](http://deeptools.readthedocs.org/en/latest/content/example_usage.html).

Please see also the [FAQ](http://deeptools.readthedocs.org/en/latest/content/help_faq.html), which we update regularly.
Our [Gallery](http://deeptools.readthedocs.org/en/latest/content/example_gallery.html) may give you some more ideas about the scope of deepTools.

For more specific **troubleshooting, feedback, and tool suggestions**, contact us via deeptools@googlegroups.com.


-------------------------------------------------------------------------------------------------------------------

<a name="installation"/></a>
Installation
---------------

deepTools are available for:

* command line usage
* integration into Galaxy servers

Details on the installation routines can be found here.

[Linux/Mac Installation](#general)

[Galaxy installation](#galaxy)


<a name="general"/></a>
### Linux/Mac Installation

The easiest way to install deepTools is by using python `pip` or `easy_install tools`:

Requirements: Python 2.7, numpy, scipy (http://www.scipy.org/install.html), bx-python, pysam, and pyBigWig

Commands:

      $ pip install deeptools --user
Done.


__Using anaconda:__

[Anaconda](https://www.continuum.io/downloads) already comes with scipy, numpy and matplotlib, making installation very quick. To install using either Anaconda or Miniconda:

    $ conda install -c bioconda deeptools

Note that deepTools does not (yet) work with python3.

__Another option is to clone the repository:__
	
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

To install into your home directory, use:

	$ python setup.py install --user

<a name="galaxy"/></a>
### Galaxy Installation

deepTools can be easily integrated into [Galaxy](http://galaxyproject.org). All wrappers and dependencies are 
available in the [Galaxy Tool Shed](http://toolshed.g2.bx.psu.edu/view/bgruening/deeptools).


#### Installation via Galaxy API (recommended)

At first generate an [API Key](http://wiki.galaxyproject.org/Admin/API#Generate_the_Admin_Account_API_Key) for your admin 
user and run the the installation script:

	python ./scripts/api/install_tool_shed_repositories.py --api YOUR_API_KEY -l http://localhost --url http://toolshed.g2.bx.psu.edu/ -o bgruening -r <revision> --name deeptools --tool-deps --repository-deps --panel-section-name deepTools

The -r argument specifies the version of deepTools. You can get the latest revsion number from the test tool shed or with the following command:

	hg identify http://toolshed.g2.bx.psu.edu/view/bgruening/deeptools

You can watch the installation status under: Top Panel → Admin → Manage installed tool shed repositories


#### Installation via web browser

- go to the [admin page](http://localhost:8080/admin)
- select *Search and browse tool sheds*
- Galaxy tool shed → Sequence Analysis → deeptools
- install deeptools

remember: for support, questions, or feature requests contact: deeptools@googlegroups.com

------------------------------------

This tool suite is developed by the [Bioinformatics Facility](http://www1.ie-freiburg.mpg.de/bioinformaticsfac) at the [Max Planck Institute for Immunobiology and Epigenetics, Freiburg](http://www1.ie-freiburg.mpg.de/).

[Documentation](http://deeptools.readthedocs.org/en/latest/index.html) | [deepTools Galaxy](http://deeptools.ie-freiburg.mpg.de) | [FAQ](http://deeptools.readthedocs.org/en/latest/content/help_faq.html)
