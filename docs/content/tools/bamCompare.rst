bamCompare
===========

``bamCompare`` can be used to generate a :ref:`bigWig` or :ref:`bedGraph` file based on **two BAM** files that are compared to each other while being simultaneously normalized for sequencing depth.

.. image:: ../../images/norm_IGVsnapshot_indFiles.png

If you are not familiar with BAM, bedGraph and bigWig formats, you can read up on that in our :doc:`../help_glossary`

The basic algorithm works proceeds in two steps:

1. Per-sample scaling / depth Normalization:

   - If scaling is used (using the SES or read counts method), appropriate scaling
     factors are determined to account for sequencing depth differences.
   - Optionally scaling can be turned off and individual samples normalized using the
     RPKM, BPM or CPM methods (or no normalization at all)

2. A per-bin calculation is performed after accounting for scaling:

   - The genome is transversed and the log2 ratio/ratio/difference/etc. for each bin of fixed width is computed.


.. argparse::
   :ref: deeptools.bamCompare.parseArguments
   :prog: bamCompare
   :nodefault:
