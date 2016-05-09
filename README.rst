======================================================================
deepTools
======================================================================
|Build Status| |Documentation Status| |PyPI version| |bioconda-badge|

User-friendly tools for exploring deep-sequencing data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

deepTools addresses the challenge of handling the large amounts of data
that are now routinely generated from DNA sequencing centers. deepTools
contains useful modules to process the mapped reads data for multiple
quality checks, creating **normalized coverage files** in standard
bedGraph and bigWig file formats, that allow comparison between
different files (for example, treatment and control). Finally, using
such normalized and standardized files, deepTools can create many
publication-ready **visualizations** to identify enrichments and for
functional annotations of the genome.

For support, questions, or feature requests contact:
deeptools@googlegroups.com

Citation:
^^^^^^^^^

Ramírez F, Ryan DP, Grüning B, Bhardwaj V, Kilpert F, Richter AS, Heyne
S, Dündar F, Manke T. `deepTools2: a next generation web server for
deep-sequencing data
analysis. <https://nar.oxfordjournals.org/content/early/2016/04/12/nar.gkw257.abstract>`__
Nucleic Acids Research. 2016 Apr 13:gkw257.

Documentation:
^^^^^^^^^^^^^^

Our `documentation <http://deeptools.readthedocs.org/>`__ contains more
details on the `individual tool scopes and
usages <http://deeptools.readthedocs.org/en/latest/content/list_of_tools.html>`__
and an `introduction to our deepTools Galaxy web
server <http://deeptools.readthedocs.org/en/latest/content/help_galaxy_intro.html>`__
including `step-by-step
protocols <http://deeptools.readthedocs.org/en/latest/content/example_usage.html>`__.

    Please see also the
    `FAQ <http://deeptools.readthedocs.org/en/latest/content/help_faq.html>`__,
    which we update regularly. Our
    `Gallery <http://deeptools.readthedocs.org/en/latest/content/example_gallery.html>`__
    may give you some more ideas about the scope of deepTools.

    For more specific **troubleshooting, feedback, and tool
    suggestions**, contact us via deeptools@googlegroups.com.

--------------

Installation
^^^^^^^^^^^^

deepTools are available for:

-  Command line usage (via pip/anaconda/github)
-  Integration into Galaxy servers (via toolshed/API/web-browser)

There are many easy ways to install deepTools. Details can be found
`here <https://deeptools.readthedocs.io/en/latest/content/installation.html>`__

**Install by cloning this repository:**

You can install any one of the deepTools branches on command line
(linux/mac) by cloning this git repository :

::

    $ git clone https://github.com/fidelram/deepTools
    $ cd deepTools
    $ python setup.py install

By default, the script will install the python library and executable
codes globally, which means you need to be root or administrator of the
machine to complete the installation. If you need to provide a
nonstandard install prefix, or any other nonstandard options, you can
provide many command line options to the install script.

::

    $ python setup.py --help

For example, to install under a specific location use:

::

    $ python setup.py install --prefix <target directory>

To install into your home directory, use:

::

    $ python setup.py install --user

**Note:** From version 2.3 onwards, deepTools support **python3**. In
case of any problems running with python3/python2, contact our user
group : deeptools@googlegroups.com.

--------------

This tool suite is developed by the `Bioinformatics
Facility <http://www1.ie-freiburg.mpg.de/bioinformaticsfac>`__ at the
`Max Planck Institute for Immunobiology and Epigenetics,
Freiburg <http://www1.ie-freiburg.mpg.de/>`__.

`Documentation <http://deeptools.readthedocs.org/en/latest/index.html>`__
\| `deepTools Galaxy <http://deeptools.ie-freiburg.mpg.de>`__ \|
`FAQ <http://deeptools.readthedocs.org/en/latest/content/help_faq.html>`__

.. |Build Status| image:: https://travis-ci.org/fidelram/deepTools.svg?branch=master
   :target: https://travis-ci.org/fidelram/deepTools
.. |Documentation Status| image:: https://readthedocs.org/projects/deeptools/badge/
   :target: http://deeptools.readthedocs.org/
.. |PyPI version| image:: https://badge.fury.io/py/deeptools.svg
   :target: https://badge.fury.io/py/deeptools
.. |bioconda-badge| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
   :target: http://bioconda.github.io
