bamPEFragmentSize
=================

.. argparse::
   :ref: deeptools.bamPEFragmentSize.parse_arguments
   :prog: bamPEFragmentSize
   :nodefault:

Example usage
^^^^^^^^^^^^^^

.. code:: bash

    $ deepTools2.0/bin/bamPEFragmentSize \
    -hist fragmentSize.png \
    -T "Fragment size of PE RNA-seq data" \
    --maxFragmentLength 1000 \
    -b testFiles/RNAseq_sample1.bam testFiles/RNAseq_sample2.bam \
    testFiles/RNAseq_sample3.bam testFiles/RNAseq_sample4.bam \
    -samplesLabel sample1 sample2 sample3 sample4


## Output

.. code:: bash

    BAM file : testFiles/RNAseq_sample1.bam

    Sample size: 10815


    Fragment lengths:
    Min.: 0.0
    1st Qu.: 311.0
    Mean: 8960.68987517
    Median: 331.0
    3rd Qu.: 362.0
    Max.: 53574842.0
    Std: 572421.46625

    Read lengths:
    Min.: 20.0
    1st Qu.: 101.0
    Mean: 99.1621821544
    Median: 101.0
    3rd Qu.: 101.0
    Max.: 101.0
    Std: 9.16567362755

    BAM file : testFiles/RNAseq_sample2.bam

    Sample size: 6771


    Fragment lengths:
    Min.: 43.0
    1st Qu.: 148.0
    Mean: 176.465071629
    Median: 164.0
    3rd Qu.: 185.0
    Max.: 500.0
    Std: 53.733877263

    ......(output truncated)


.. image:: ../../images/test_plots/ExampleFragmentSize.png
