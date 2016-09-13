multiBamSummary
================

.. contents:: 
    :local:

.. argparse::
   :ref: deeptools.multiBamSummary.parse_arguments
   :prog: multiBamSummary
   :nodefault:


Example
^^^^^^^

The default output of ``multiBamSummary`` (a compressed ``numpy`` array: `*.npz`) can be visualized using :doc:`plotCorrelation` or :doc:`plotPCA`.

The optional output (``--outRawCounts``) is a simple tab-delimited file that can be used with any other program. The first three columns define the region of the genome for which the reads were summarized.

.. code:: bash

    $ deepTools2.0/bin/multiBamSummary bins \
     --bamfiles testFiles/*bam \ # using all BAM files in the folder
     --minMappingQuality 30 \
     --region 19 \ # limiting the binning of the genome to chromosome 19
     --labels H3K27me3 H3K4me1 H3K4me3 HeK9me3 input \
     -out readCounts.npz --outRawCounts readCounts.tab

     $ head readCounts.tab 
     #'chr'	'start'	'end'	'H3K27me3'	'H3K4me1'	'H3K4me3'	'HeK9me3'	'input'
     19	10000	20000	0.0	0.0	0.0	0.0	0.0
     19	20000	30000	0.0	0.0	0.0	0.0	0.0
     19	30000	40000	0.0	0.0	0.0	0.0	0.0
     19	40000	50000	0.0	0.0	0.0	0.0	0.0
     19	50000	60000	0.0	0.0	0.0	0.0	0.0
     19	60000	70000	1.0	1.0	0.0	0.0	1.0
     19	70000	80000	0.0	1.0	7.0	0.0	1.0
     19	80000	90000	15.0	0.0	0.0	6.0	4.0
     19	90000	100000	73.0	7.0	4.0	16.0	5.0
