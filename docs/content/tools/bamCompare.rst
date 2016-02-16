bamCompare
===========

``bamCompare`` can be used to generate a :ref:`bigWig` or :ref:`bedGraph` file based on **two BAM** files that are compared to each other while being simultaneously normalized for sequencing depth.

.. image:: ../../images/norm_IGVsnapshot_indFiles.png

If you are not familiar with BAM, bedGraph and bigWig formats, you can read up on that in our :doc:`../help_glossary`

.. argparse::
   :ref: deeptools.bamCompare.parseArguments
   :prog: bamCompare

.. warning:: The filtering by deepTools is done *after* the scaling factors are calculated!

.. warning:: If you know that your files will be strongly affected by the kind of filtering you would like to apply (e.g., removal of duplicates with ``--ignoreDuplicates`` or ignoring reads of low quality) then consider removing those reads *beforehand*. 
