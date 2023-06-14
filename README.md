# deepTools
[![Documentation Status](https://readthedocs.org/projects/deeptools/badge/)](http://deeptools.readthedocs.org/) 
[![PyPI Version](https://img.shields.io/pypi/v/deeptools.svg?style=plastic)](https://pypi.org/project/deepTools/) 
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/deeptools/README.html)
[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=deeptools_compute_matrix)
![test](https://github.com/deeptools/deepTools/actions/workflows/test.yml/badge.svg)
![planemo](https://github.com/deeptools/deepTools/actions/workflows/planemo.yml/badge.svg)


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
