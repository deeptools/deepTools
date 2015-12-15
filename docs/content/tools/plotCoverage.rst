plotCoverage
===========

.. argparse::
   :ref: deeptools.plotCoverage.parse_arguments
   :prog: plotCoverage

   
Usage Example:
~~~~~~~~~~~~~~

plotCoverage -b H3K4Me1.bam H3K4Me3.bam H3K27Me3.bam H3K9Me3.bam --plotFile example_coverage -n 1000000 -p 5 --plotTitle "example_coverage" --ignoreDuplicates --minMappingQuality 10 --region 19
