bamPEFragmentSize
=================

.. argparse::
   :ref: deeptools.bamPEFragmentSize.parse_arguments
   :prog: bamPEFragmentSize
   
.. note:: This tool accepts only one BAM file at a time.

Example usage
^^^^^^^^^^^^^^

.. code:: bash

    $ deepTools2.0/bin/bamPEFragmentSize \
    -hist fragmentSize.png \
    -T "Fragment size of PE RNA-seq data" \
    testFiles/RNAseq.bam 
 
     Sample size: 12850

    Fragment lengths:
        Min.: 0.0
        1st Qu.: 313.0
        Mean: 2597.03237354
        Median: 357.0
        3rd Qu.: 2726.0
        Max.: 384622.0
        Std: 7066.11863701

    Read lengths:
        Min.: 20.0
        1st Qu.: 101.0
        Mean: 99.4182101167
        Median: 101.0
        3rd Qu.: 101.0
        Max.: 101.0
        Std: 7.64455778462

.. image:: ../../images/test_plots/ExampleFragmentSize.png