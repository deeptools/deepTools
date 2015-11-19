deepTools
=========

.. image:: ../../images/start_collage.png

For support, questions, or feature requests contact:
deeptools@googlegroups.com

The main reason why deepTools was started is the simple fact that in
2011 we could not find tools that met all our needs for NGS data
analysis. While there were individual tools for separate tasks, we
wanted software that would fulfill *all* of the following criteria:

-  **efficiently extract reads from BAM files** and perform various
   computations on them
-  **turn BAM files of aligned reads into bigWig files** using different
   normalization strategies
-  make use of **multiple processors** (speed!)
-  generation of **highly customizable images** (change colours, size,
   labels, file format etc.)
-  enable **customized down-stream analyses** which requires that every
   data set that is being produced can be stored by the user
-  **modular approach** - compatibility, flexibility, scalability (i.e.
   we can add more and more modules making use of established methods)

The flow chart below depicts the different tool modules that are
currently available within deepTools (deepTools modules are written in
bold red and black font).

.. image:: ../../images/start_workflow.png

If  the file names in the figure mean nothing to you,
please make sure to check our `Glossary`_.

+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| tool          | type          | input files                       | main output file(s)                    | application                                                                  |
+===============+===============+===================================+========================================+==============================================================================+
| bamCorrelate  | QC            | 2 or more BAM                     | clustered heatmap                      | Pearson or Spearman correlation between read distributions                   |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| bamFingerprin | QC            | 2 BAM                             | 1 diagnostic plot                      | assess enrichment strength of a ChIP sample                                  |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| computeGCbias | QC            | 1 BAM                             | 2 diagnostic plots                     | calculate the exp. and obs. GC distribution of reads                         |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| correctGCbias | QC            | 1 BAM, output from computeGCbias  | 1 GC-corrected BAM                     | obtain a BAM file with reads distributed according to the genome’s GC conten |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| bamCoverage   | normalization | BAM                               | bedGraph or bigWig                     | obtain the normalized read coverage of a single BAM file                     |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| bamCompare    | normalization | 2 BAM                             | bedGraph or bigWig                     | normalize 2 BAM files to each other (e.g. log2ratio, difference)             |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| computeMatrix | visualization | 1 bigWig, 1 BED                   | zipped file for heatmapper or profiler | compute the values needed for heatmaps and summary plots                     |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| heatmapper    | visualization | computeMatrix output              | heatmap of read coverages              | visualize the read coverages for genomic regions                             |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+
| profiler      | visualization | computeMatrix output              | summary plot (“meta-profile”)          | visualize the average read coverages over a group of genomic regions         |
+---------------+---------------+-----------------------------------+----------------------------------------+------------------------------------------------------------------------------+


Contents:
---------
.. toctree::
   :glob:
   :maxdepth: 2

   content/list_of_tools
   *
   Help/*
   technical_documentation/*
   tools_details/*
   example_workflows/*
   api_tutorial


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

About
-----

Please cite deepTools as follows:
	  Fidel Ramírez, Friederike Dündar, Sarah Diehl, Björn A. Grüning, and Thomas Manke.
	  `deepTools: a flexible platform for exploring deep-sequencing data. <http://nar.oxfordjournals.org/content/early/2014/05/05/nar.gku365.abstract>`_
	  Nucl. Acids Res., 2014
	  doi:10.1093/nar/gku365

.. image:: ../../images/logo_mpi-ie.jpg
	  
This tool suite is developed by the `Bioinformatics Facility <http://www1.ie-freiburg.mpg.de/bioinformaticsfac>`_ at the
`Max Planck Institute for Immunobiology and Epigenetics,
Freiburg <http://www1.ie-freiburg.mpg.de/>`_.

`deepTools Galaxy`_ | `get in touch`_!

.. _Glossary: <https://github.com/fidelram/deepTools/wiki/Glossary>
.. _deepTools Galaxy: <http://deeptools.ie-freiburg.mpg.de>
.. _get in touch: <deeptools@googlegroups.com>