Example usage
=============

.. toctree::
   :maxdepth: 1

   step_by_step_protocols
   gallery
   api_tutorial

How we use deepTools 
--------------------

deepTools started off as a tool package for ChIP-seq analyses, which is
why you find many ChIP-seq examples in our documentation.
`Here <https://https.google.com/file/d/0B8DPnFM4SLr2UjdYNkQ0dElEMm8/edit?usp=sharing>`__
are slides that we used for teaching at the University of Freiburg with
more details on the deepTools usage and aims in regard to ChIP-seq.
However, while some tools, such as `bamFingerprint`, specifically
address ChIP-seq-issues, the majority of tools is widely applicable
to NGS data, including RNA-seq.

So, how does a typical, basic ChIP-seq workflow look like for us?

.. image:: ../images/start_workflow.png

As shown in the flow chart above, our work usually begins with one or
more `FASTQ <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-fastq>`__
file(s) of deeply-sequenced samples. After a first quality control using
`FASTQC <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`__,
we align the reads to the reference genome, e.g. using
`bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`__.
We then use deepTools to assess the quality of the aligned reads:

#. **Correlation between BAM files** (``bamCorrelate``). This is a very basic test to see whether
   the sequenced and aligned reads meet your expectations. We use this
   check to assess the reproducibility - either between replicates
   and/or between different experiments that might have used the same
   antibody or the same cell type etc. For instance, replicates should
   correlate better than differently treated samples.
#. **GC bias check** (``computeGCbias``). Many sequencing protocols
   require several rounds of PCR-based amplification of the DNA to be
   sequenced. Unfortunately, most DNA polymerases used for PCR introduce
   significant GC biases as they prefer to amplify GC-rich templates.
   Depending on the sample (preparation), the GC bias can vary
   significantly and we routinely check its extent. In case we need to
   compare files with different GC biases, we use the *correctGCbias*
   module to match the GC bias.
   See the paper by `Benjamini and
   Speed <http://nar.oxfordjournals.org/content/40/10/e72>`__ for many
   insights into this problem.
#. **Assessing the ChIP strength**. We do this quality control to get a
   feeling for the signal-to-noise ratio in samples from ChIP-seq
   experiments. It is based on the insights published by `Diaz et
   al. <http://www.degruyter.com/view/j/sagmb.2012.11.issue-3/1544-6115.1750/1544-6115.1750.xml>`__.

Once we're satisfied by the basic quality checks, we normally **convert**
the large
`BAM <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bam>`__
files into a leaner data format, typically
`bigWig <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bigwig>`__.
bigWig files have several advantages over BAM files that mainly stem
from their significantly decreased size:

-  useful for data sharing and storage
-  intuitive visualization in Genome Browsers (e.g.
   `IGV <http://www.broadinstitute.org/igv/>`__)
-  more efficient downstream analyses are possible

The deepTools modules ``bamCompare`` and ``bamCoverage`` do not only allow
the simple conversion from BAM to bigWig (or
bedGraph for that matter), **the main reason why we developed those tools was
that we wanted to be able to *normalize* the read coverages** so that we
could compare different samples despite differences in sequencing depth,
GC biases and so on.

Finally, once all the files have passed our visual inspections, the fun
of downstream analyses with ``computeMatrix``, ``heatmapper`` and ``profiler`` can begin!