*********
deepTools
*********

This suite of tools provides useful routines to deal with mapped reads.
This includes converting bams to bigwig or bedgraph coverage files, run a 
GC correction and plot heatmaps and profiles.

======================
Install from source
===================

The easiest way to install deppTools is by downloading the
source file and using python distutils tools::

 $ wget https://github.com/fidelram/deepTools/blob/master/dist/deepTools-1.5.tar.gz
 $ tar -xzvf deepTools-1.5.tar.gz

Then go to the unpacked directory and
and simply run the install script::

 $ cd deepTools-1.5
 $ python setup.py install

By default, the script will install python library and executable
codes globally, which means you need to be root or administrator of
the machine to complete the installation. If you need to
provide a nonstandard install prefix, or any other nonstandard
options, you can provide many command line options to the install
script. Use the

 $ python setup.py --help

To install under a specific location use::

 $ python setup.py install --prefix <target directory>

