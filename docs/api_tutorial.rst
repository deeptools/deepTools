.. _api:

deepTools API example
=====================

deepTools consists of several command line and galaxy wrappers for summarising
the information on Next Generation Sequencing data that can be mapped to a reference
genome. However, the engine powering the deepTools commands can be used through the API for other
purposes as well. The following is a short overview of the most useful methods and classes
from deepTools.


Finding the coverage over a region
----------------------------------

With deepTools the read coverage over multiple genomic regions and
multiple files can be computed quite quickly using multiprocessing. First, we
start with simple example that is later extended to cover the use of multiprocessing.
In this example we compute the coverage of reads over a small region for bins of 50bp. For
this we need the :class:`deeptools.countReadsPerBin` class.


.. code:: python

    import deeptools.countReadsPerBin


We also need a bam file containing the aligned reads to a reference genome. The bam file must
be indexed to allow quick access to the reads falling into the regions of interest.

.. code:: python

    bam_file = "file.bam"

Now, the countReadsPerBin object can be initialized.
The first argument of the constructor is a list of bam files, which in this case is
just one. We are going to use a
binLength of 50 bp that is computed every 50 bp pairs, in other words, in
consecutive bins of 50bp. Overlapping bin
coverages can be used by setting a `stepSize` smaller than the bin length.

.. code:: python

    cr = countReadsPerBin.CountReadsPerBin([bam_file], binLength=50, stepSize=50)


Now, we can compute the coverage over a region in chromosome 2 from position 0
to 1000.

.. code:: python

    cr.count_reads_in_region('chr2L', 0, 1000)

.. parsed-literal::

    array([[ 2.],
           [ 3.],
           [ 1.],
           [ 2.],
           [ 3.],
           [ 2.],
           [ 4.],
           [ 3.],
           [ 2.],
           [ 3.],
           [ 4.],
           [ 6.],
           [ 4.],
           [ 2.],
           [ 2.],
           [ 1.]])

The result is a numpy array in which each row corresponds to a bin and each column corresponds
to a bam file. In this case only one bam file was used.

Filtering regions
-----------------

If reads want to be filtered the relevant options simply
need to be passed to the constructor. In the following code, the reads are filtered
using a minimum mapping quality of 20, and by including only reads that map to the forward
strand using a sam flag (samFlag_exclude=16, where 16 is the value for reverse reads).
Furthermore, duplicated reads are ignored.

.. code:: python

    cr = countReadsPerBin.CountReadsPerBin([bam_file], binLength=50, stepSize=50,
                                            minMappingQuality=20,
                                            samFlag_exclude=16,
                                            ignoreDuplicates=True
                                            )
    cr.count_reads_in_region('chr2L', 1000000, 1001000)

.. parsed-literal::

    array([[ 1.],
           [ 1.],
           [ 0.],
           [ 0.],
           [ 0.],
           [ 0.],
           [ 2.],
           [ 3.],
           [ 1.],
           [ 0.],
           [ 1.],
           [ 2.],
           [ 0.],
           [ 0.],
           [ 1.],
           [ 2.],
           [ 1.],
           [ 0.],
           [ 0.],
           [ 0.]])

Sampling the genome
-------------------

Instead of consecutive bins as in the previous cases, a genome can
simply be sampled. This is useful to estimate some values,
like depth of sequencing, without having to look at the complete genome. In the following example,
10.000 positions of size 1 bp are going to be queried from three bam files to compute the average depth of sequencing.
For this we set the numberOfSamples parameter in the object constructor. The `skipZeros` parameter
is added such that regions that in all bam files do not have any reads are excluded. Usually, those
regions are repetive which are often excluded from the read mapping. The `run()` method is
used instead of `count_reads_in_region`.

.. code:: python

    cr = countReadsPerBin.CountReadsPerBin([bam_file1, bam_file2, bam_file3],
                                            binLength=1, numberOfSamples=10000,
                                            numberOfProcessors=10,
                                            skipZeros=True)
    sequencing_depth = cr.run()
    print sequencing_depth.mean(axis=0)

.. parsed-literal::
    [  1.98923924   2.43743744  22.90102603]


The `run()` method splits the computation of the coverage over 10 processors and aggregates
the results. When the parameter number of samples is used the regions selected
for the computation of the coverage are not random. Instead, the genome is split into 'number-of-samples'
equal parts and at the start of each part is then queried for the coverage. If truly random values are
required is recommended to pass a bed file to the constructor containing the regions to be sampled.


Now it is possible to make some diagnostic plots from the results:

.. code:: python

    fig, axs = plt.subplots(1, 2, figsize=(15,5))
    # plot coverage
    for col in res.T:
        axs[0].plot(np.bincount(col.astype(int)).astype(float)/total_sites)
        csum = np.bincount(col.astype(int))[::-1].cumsum()
        axs[1].plot(csum.astype(float)[::-1] / csum.max())
    axs[0].set_xlabel('coverage')
    axs[0].set_ylabel('fraction of bases sampled')
    # plot cumulative coverage

    axs[1].set_xlabel('coverage')
    axs[1].set_ylabel('fraction of bases sampled >= coverage')


.. image:: images/plot_coverage.png


Computing the FRiP score
------------------------

The FRiP score is defined as the fraction of reads that fall into a peak and is 
often used as a measure of ChIP-seq quality. For this example we
need a  bed file containing the peak regions. Such file is
usually computed using a peak caller. Also, two bam files are
going to be used that correspond to two biological replicates.

.. code:: python

    bed_file = open("peaks.bed", 'r')
    cr = countReadsPerBin.CountReadsPerBin([bam_file1, bam_file2],
                                            bedFile=bed_file,
                                            numberOfProcessors=10)
    reads_at_peaks = cr.run()
    print reads_at_peaks

.. parsed-literal::

    array([[ 322.,  248.],
           [ 231.,  182.],
           [ 112.,  422.],
           ..., 
           [ 120.,   76.],
           [ 235.,  341.],
           [ 246.,  265.]])


The result is a numpy array having as rows each of the peak regions and as columns each of the bam files.

.. code:: python

    reads_at_peaks.shape


.. parsed-literal::

    (6295, 2)

Now, the total number of reads per falling within the peaks, per bam file, is computed:

.. code:: python

    total = reads_at_peaks.sum(axis=0)

Next, we need to find the total number of mapped reads in each of the bam files. For
this we use the pysam module.

.. code:: python

    import pysam
    bam1 = pysam.AlignmentFile(bam_file1)
    bam2 = pysam.AlignmentFile(bam_file2)

Now, `bam1.mapped` and `bam2.mapped` contain the total number of mapped
reads in each of the bam files respectively.

Finally, we can compute the FRiP score:

.. code:: python

    frip1 = float(total[0]) / bam1.mapped
    frip2 = float(total[1]) / bam2.mapped
    print frip1, frip2

.. parsed-literal::

    0.170030741997, 0.216740390353



Using mapReduce to sample paired-end fragment lengths
------------------------------------------------------

deepTools internally uses a map-reduce strategy in which a computation is split into smaller
parts that are sent to different processors which is subsequently integrated. The following
example is based on the code available for `bamPEFragmentSize.py`

In this case retrieve the reads from a bam file and collect the
fragment length. Reads are retrieved using pysam, and the `read` object returned
contains the `template_length` attribute which is the number of bases from the
leftmost mapped base to the rightmost mapped base in the read pair.

First, we will create a function that can collect fragment lengths over a genomic
position from a bam file. Because later we will call this function using
mapReduce the function accepts only one argument that is
a tuple in which the first three parameters are set to
chromosome name, start and end. The next parameter is the bam file name.

.. code:: python

    import pysam
    import numpy as np
    def get_fragment_length(args):
        chrom, start, end, bam_file_name = args
        bam = pysam.Aligmementfile(bam_file_name)
        f_lens_list = []
        for fetch_start in range(start, end, 1e6):
            # simply get the reads over a region of 10000 bp
            fetch_end = min(end, start + 10000)

            f_lens_list.append(np.array([abs(read.template_length)
                                  for read in bam.fetch(chrom, fetch_start, fetch_end)
                                  if read.is_proper_pair and read.is_read1]))

        # concatenate all results
        return np.concatenate(fragment_lengths)


Now, we can use `mapReduce` to call this function and compute fragment lengths
over the whole genome. mapReduce needs to know the chromosome sizes which
can be easily retrieved from the bam file. Furthermore, it needs to know
the size of the region that is sent to each processor. For this
example, a region of 10 million bp is sent to each processor where
the function just defined (get_fragment_length) is going to be called. In other
words, each processor executes the same get_fragment_length function to collect data over
a 10 million bp. The arguments to mapReduce are the list of arguments sent to the function, besides
the first obligatory three (chrom start, end). In this case only one extra argument is passed
to the function, the bam file name. The next two positional arguments are the name of the function to call
(`get_fragment_length`) the the chromosome sizes.

.. code:: python

    import deeptools.mapReduce
    bam = pysam.Aligmentfile(bamFile)
    chroms_sizes = zip(bam.references, bam.lengths)

    result = mapReduce.mapReduce((bam_file_name, ),
                                  get_fragment_length
                                  chrom_sizes,
                                  genomeChunkLength=10000000,
                                  numberOfProcessors=20,
                                  verbose=True)

    fragment_lengths =  np.concatenate(result)

    print "mean fragment length {}".format(fragment_lengths.mean()"
    print "median fragment length {}".format(np.median(fragment_lengths)"


.. parsed-literal::

    0.170030741997, 0.216740390353
