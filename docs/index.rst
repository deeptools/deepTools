=====================================================
deepTools: *tools for exploring deep sequencing data*
=====================================================

.. image:: images/start_collage.png

deepTools is a **suite of python tools** particularly developed for the 
efficient analysis of high-throughput sequencing data, such as ChIP-seq, RNA-seq or
MNase-seq.

There are 2 ways for using deepTools:

* **command line usage** -- simply download and install the tools
* **Galaxy usage** --  if you would like to see our tools at work, you can test them on our public `deepTools Galaxy server <http://deeptools.ie-freiburg.mpg.de>`_.

The flow chart below depicts the different tool modules that are
currently available (deepTools modules are written in
bold red and black font).

.. image:: images/start_workflow.png

If the file names in the figure mean nothing to you,
please make sure to check our :doc:`content/help_glossary`.


Contents:
---------
.. toctree::
   :maxdepth: 1

   content/installation
   content/list_of_tools
   content/example_usage
   content/api
   content/changelog
   content/help_galaxy_intro
   content/help_faq
   content/help_faq_galaxy
   content/help_glossary


While developing deepTools, we continuously strive to create software
that fulfills the following criteria:

-  **efficiently extract reads from BAM files** and perform various
   computations on them
-  **turn BAM files of aligned reads into bigWig files** using different
   normalization strategies
-  make use of **multiple processors** (speed!)
-  generation of **highly customizable images** (change colours, size,
   labels, file format, etc.)
-  enable **customized down-stream analyses**, meaning that every
   data set created can be stored by the user
-  **modular approach** - compatibility, flexibility, scalability (i.e.
   we can add more and more modules and make use of established methods)

For support, questions, or feature requests contact:
deeptools@googlegroups.com

About
-----

Please cite deepTools as follows:
	  Fidel Ramírez, Friederike Dündar, Sarah Diehl, Björn A. Grüning, and Thomas Manke.
	  `deepTools: a flexible platform for exploring deep-sequencing data. <http://nar.oxfordjournals.org/content/early/2014/05/05/nar.gku365.abstract>`_
	  Nucl. Acids Res., 2014
	  doi:10.1093/nar/gku365
	  
.. image:: images/logo_mpi-ie.jpg
	  
This tool suite is developed by the `Bioinformatics Facility <http://www1.ie-freiburg.mpg.de/bioinformaticsfac>`_ at the
`Max Planck Institute for Immunobiology and Epigenetics,
Freiburg <http://www1.ie-freiburg.mpg.de/>`_.

