# deepTools
[![Build Status](https://travis-ci.org/deeptools/deepTools.svg?branch=master)](https://travis-ci.org/deeptools/deepTools) [![Documentation Status](https://readthedocs.org/projects/deeptools/badge/)](http://deeptools.readthedocs.org/) [![PyPI Version](https://img.shields.io/pypi/v/deeptools.svg?style=plastic)](https://pypi.org/project/deepTools/) 
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/deeptools/README.html)

## User-friendly tools for exploring deep-sequencing data

deepTools addresses the challenge of handling the large amounts of data that are now routinely generated from DNA sequencing centers. deepTools contains useful modules to process the mapped reads data for multiple quality checks, creating **normalized coverage files** in standard bedGraph and bigWig file formats, that allow comparison between different files (for example, treatment and control). Finally, using such normalized and standardized files, deepTools can create many publication-ready  **visualizations** to identify enrichments and for functional annotations of the genome.

For support or questions please post to [Biostars](http://biostars.org). For bug reports and feature requests please open an issue [on github](http://github.com/deeptools/deeptools).


### Citation:

Ramírez F, Ryan DP, Grüning B, Bhardwaj V, Kilpert F, Richter AS, Heyne S, Dündar F, Manke T. [deepTools2: a next generation web server for deep-sequencing data analysis.](https://nar.oxfordjournals.org/content/early/2016/04/12/nar.gkw257.abstract) Nucleic Acids Research. 2016 Apr 13:gkw257.

### Documentation:

Our [documentation](http://deeptools.readthedocs.org/) contains more details on the [individual tool scopes and usages](http://deeptools.readthedocs.org/en/latest/content/list_of_tools.html) and an [introduction to our deepTools Galaxy web server](http://deeptools.readthedocs.org/en/latest/content/help_galaxy_intro.html) including [step-by-step protocols](http://deeptools.readthedocs.org/en/latest/content/example_usage.html).

>Please see also the [FAQ](http://deeptools.readthedocs.org/en/latest/content/help_faq.html), which we update regularly.
Our [Gallery](http://deeptools.readthedocs.org/en/latest/content/example_gallery.html) may give you some more ideas about the scope of deepTools.

>For more specific **troubleshooting, feedback, and tool suggestions**, please post [to Biostars](http://biostars.org).


-------------------------------------------------------------------------------------------------------------------

### Installation

deepTools are available for:

* Command line usage (via pip/anaconda/github)
* Integration into Galaxy servers (via toolshed/API/web-browser)

There are many easy ways to install deepTools. Details can be found [here](https://deeptools.readthedocs.io/en/latest/content/installation.html)

**Install by cloning this repository:**

You can install any one of the deepTools branches on command line (linux/mac) by cloning this git repository :

	$ git clone https://github.com/deeptools/deepTools
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

deepTools can be easily integrated into [Galaxy](http://galaxyproject.org). Please see the [installation instructions in our documentation](http://deeptools.readthedocs.io/en/latest/content/installation.html#galaxy-installation) for further details.

**Note:** From version 2.3 onwards, deepTools support **python3**.

------------------------------------

This tool suite is developed by the [Bioinformatics Facility](http://www1.ie-freiburg.mpg.de/bioinformaticsfac) at the [Max Planck Institute for Immunobiology and Epigenetics, Freiburg](http://www1.ie-freiburg.mpg.de/).

[Documentation](http://deeptools.readthedocs.org/en/latest/index.html) | [deepTools Galaxy](http://deeptools.ie-freiburg.mpg.de) | [FAQ](http://deeptools.readthedocs.org/en/latest/content/help_faq.html)
