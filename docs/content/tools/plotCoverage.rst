plotCoverage
============

.. argparse::
   :ref: deeptools.plotCoverage.parse_arguments
   :prog: plotCoverage

mean value indicated in the left plot as a proxy for seq. depth
what is the fraction of the genome that has a depth of sequencing of 2 

.. image:: ../../images/plotCoverage_annotated.png

Usage Example:
~~~~~~~~~~~~~~

.. code:: bash
	
   $ plotCoverage -b H3K4Me1.bam H3K4Me3.bam H3K27Me3.bam H3K9Me3.bam
      --plotFile example_coverage
      -n 1000000
      --plotTitle "example_coverage" \ 
      --ignoreDuplicates \
      --minMappingQuality 10\ 
      --region 19
      
Example plot
~~~~~~~~~~~~~~~~~~~~

.. image:: ../../images/test_plots/ExamplePlotCoverage.png

As you can see, the coverage of our test data sets is very poor -- on average, there is fewer than 1 read per bp! 