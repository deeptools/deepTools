computeMatrixOperations
=======================

.. contents:: 
    :local:

.. argparse::
   :ref: deeptools.computeMatrixOperations.parse_arguments
   :prog: computeMatrixOperations
   :nodefault:

Details
^^^^^^^

``computeMatrixOperations`` can perform a variety of operations on one or more files produced by ``computeMatrix`` (N.B., the output is always written to a new file):

+----------------+--------------------------------------------------------------------------------------------------------------------------+
+ **Subcommand** | **What it does**                                                                                                         |
+----------------+--------------------------------------------------------------------------------------------------------------------------+
+ info           | Prints out the sample and region group names in the order in which they appear.                                          |
+----------------+--------------------------------------------------------------------------------------------------------------------------+
+ subset         | Subsets a file by the desired samples/region group names. This can also change the order of these samples/region groups. |
+----------------+--------------------------------------------------------------------------------------------------------------------------+
+ filterStrand   | Filters the file to only include regions annotated as being on a particular strand.                                      |
+----------------+--------------------------------------------------------------------------------------------------------------------------+
+ rbind          | Concatenates multiple matrices together, top to bottom.                                                                  |
+----------------+--------------------------------------------------------------------------------------------------------------------------+
+ cbind          | Merges multiple matrices, left to right.                                                                                 |
+----------------+--------------------------------------------------------------------------------------------------------------------------+
+ sort           | Sorts the given file so regions are in the order of occurence in the input BED/GTF file(s).                              |
+----------------+--------------------------------------------------------------------------------------------------------------------------+


These operations are useful when you want to run computeMatrix on multiple files (thereby keeping all of the values together) and later exclude regions/samples or add new ones. Another common use would be if you require the output of computeMatrix to be sorted to match the order of regions in the input file.

.. attention::
   As of version 3.0, computeMatrix (and therefore also computeMatrixOperations) produces output with labels present for each sample. If you run any operations on matrices output by older versions then they will be modified to be comformant with the new output, which is not backward compatible!

Examples
^^^^^^^^

Suppose that we have a strand-specific RNAseq dataset and would like to plot only the strand-specific signal across spliced transcripts. The general steps would be as follows:

1. Run `bamCoverage` on each sample twice, once with `--filterRNAstrand forward` and again with `--filterRNAstrand reverse`. This will result in twice the number of bigWig files as samples.
2. Run `computeMatrix scale-regions` with all of these bigWig files, including the `--metagene` option and a BED12 and/or a GTF file. This produces a file containing the signal separated by strand for each transcript.
3. Get the list of sample names stored in the matrix file:

.. code:: bash

    $ computeMatrixOperations info -m foo.mat.gz
    Groups:
        genes
    Samples:
        SRR648667.forward
        SRR648668.forward
        SRR648669.forward
        SRR648670.forward
        SRR648667.reverse
        SRR648668.reverse
        SRR648669.reverse
        SRR648670.reverse

4. Create two new files, each containing only the samples containing signal from a given strand.

.. code:: bash

    $ computeMatrixOperations subset -m foo.mat.gz -o forward.mat.gz --samples SRR648667.forward SRR648668.forward SRR648669.forward SRR648670.forward
    $ computeMatrixOperations subset -m foo.mat.gz -o reverse.mat.gz --samples SRR648667.reverse SRR648668.reverse SRR648669.reverse SRR648670.reverse

5. These files can then be subset to contain only transcripts on a particular strand. Note that it's best to double check that the ``--strand -`` setting produces the intended results. There are many peculiar variants of RNAseq library preparation and the settings for one type may not be appropriate for another (to check this, use the different ``--strand`` options on the same matrix and then run ``plotHeatmap``, one of them will be obviously correct and the other largely blank).

.. code:: bash

    $ computeMatrixOperations filterStrand -m forward.mat.gz -o forward.subset.mat.gz --strand -
    $ computeMatrixOperations filterStrand -m reverse.mat.gz -o reverse.subset.mat.gz --strand +

6. Finally, the files can be merged back together, head to tail. The samples are already in the correct order, as indicated by step 3.

.. code:: bash

    $ computeMatrixOperations rbind -m forward.subset.mat.gz reverse.subset.mat.gz -o merged.mat.gz

7. If desired, the transcripts can then be resorted to match the order of the input GTF file.

.. code:: bash

    $ computeMatrixOperations sort -m merged.mat.gz -o sorted.mat.gz -R genes.gtf

The resulting file can then be used with ``plotHeatmap`` or ``plotProfile``. Note that we could have skipped the subset step and run ``computeMatrix`` independently on the forward and reverse bigWig files.
