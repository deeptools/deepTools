======================================================================
deepTools: suite of common processing routies for deep sequencing data 
======================================================================

This suite of tools provides useful routines to deal with mapped reads.
This includes converting bams to bigwig or bedgraph coverage files, run a 
GC correction and plot heatmaps and profiles.

===================
Install from source
===================

The easiest way to install deepTools is by downloading the
source file and using python pip or easy_install tools:

	$ pip install deepTools

Other option is to clone the repository
	
	$ git clone https://github.com/fidelram/deepTools

Then go to the deepTools directory and
and simply run the install script:

	$ cd deepTools
	$ python setup.py install

By default, the script will install python library and executable
codes globally, which means you need to be root or administrator of
the machine to complete the installation. If you need to
provide a nonstandard install prefix, or any other nonstandard
options, you can provide many command line options to the install
script.

	$ python setup.py --help

To install under a specific location use::

	$ python setup.py install --prefix <target directory>

