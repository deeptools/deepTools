`General deepTools-related questions <#general>`__
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  `How does deepTools handle data from paired-end sequencing? <#PE>`__
-  `I just want to try out a tool, how can I optimize the computation
   time? <#compTime>`__
-  `Can I specify more than one chromosome in the --regions
   Option? <#regionsLimit>`__
-  `When should I exclude regions from computeGCbias? <#excludeGC>`__
-  `When should I use bamCoverage, when
   bamCompare? <#bamCompare_vs_bamCoverage>`__

computeMatrix/heatmapper
                        

-  `How does computeMatrix handle overlapping genome
   regions? <#compMatrix_overlap>`__
-  `Why does the maximum value in the heatmap not equal the maximum
   value in the matrix? <#hm_scaling>`__
-  `How can I increase the resolution of the heatmap? <#hmresolution>`__
-  `How can I change the automatic labels of the clusters in a k-means
   clustered heatmap? <#hmlabels>`__
-  `How can I manually specify several groups of regions (instead of
   clustering)? <#hmgroups>`__

External data
             

-  `What do I have to pay attention to when working with a draft version
   of a genome? <#draftGenome>`__
-  `How do I calculate the effective genome size for an organism that's
   not in your list? <#effGenomeSize>`__
-  `Where can I download the 2bit genome files required for
   computeGCbias? <#2bit>`__

--------------

How does deepTools handle data from paired-end sequencing?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| Generally, all the modules working with
`BAM <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bam>`__
files (*bamCorrelate, bamCoverage, bamCompare, bamFingerprint,
computeGCbias*) recognize paired-end sequencing data. You can enforce to
ignore the fragment length based on the mate pairs using the option
*doNotExtendPairedEnds* ("advanced options" in Galaxy).
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_
| 

How can I test a tool with little computation time?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When you're playing around with the tools to see what kinds of results
they will produce, you can limit the operation to one chromosome or a
specific region save computation time. In Galaxy, you will find this
under "advanced output options" → "Region of the genome to limit the
operation to"; the command line option is called "--region"
(CHR:START:END).

The following tools currently have this option:

-  `bamCorrelate <https://github.com/fidelram/deepTools/wiki/QC>`__
-  `bamFingerprint <https://github.com/fidelram/deepTools/wiki/QC>`__
-  `computeGCbias <https://github.com/fidelram/deepTools/wiki/QC>`__,
   `correctGCbias <https://github.com/fidelram/deepTools/wiki/Normalizations>`__
-  `bamCoverage <https://github.com/fidelram/deepTools/wiki/Normalizations>`__,
   `bamCompare <https://github.com/fidelram/deepTools/wiki/Normalizations>`__

It works as follows: first, the *entire* genome represented in the
`BAM <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bam>`__
file will be regarded and sampled, *then* all the regions or sampled
bins that do not overlap with the region indicated by the user will be
discarded.

| Be aware that you can limit the operation to only *one* chromosome (or
*one* specific locus on a chromosome), if you would like to limit the
operation to more than one region, see the `answer
below <#regionsLimit>`__.
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_
| 

Can I specify more than one chromosome in the --regions Option?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several programs have the option to limit the operation to a specific
region (i.e. to reduce computation time). The input must be given in the
format of chr:start:end, for example "chr10" or "chr10:456700:891000".
It is not possible to indicate more than one region, e.g. chr10, chr11 -
**this will not work**!

Here are some ideas for workarounds if you do find yourself in a
situation where you would in fact like to limit the operation to more
than one region:

-  **general workaround**: since all the tools that have the --region
   option work on BAM files, you could *filter your reads* prior to
   running the program, e.g. using intersectBed with --abam or samtools
   view. Use the resulting (smaller) BAM file for the deepTools program
   of your choice.

::

    samtools view -b -L regionsOfInterest.bed Reads.bam > ReadsOverlappingWithRegionsOfInterest.bam

::

    intersectBed -abam Reads.bam -b regionsOfInterest.bed > ReadsOverlappingWithRegionsOfInterest.bam

`computeGCbias <https://github.com/fidelram/deepTools/wiki/QC>`__ and
`bamCorrelate <https://github.com/fidelram/deepTools/wiki/QC>`__ offer
in-build solutions:

-  **bamCorrelate**: bamCorrelate has two modes, bins and BED. If you
   make use of the BED mode (for details, see
   `here <https://github.com/fidelram/deepTools/wiki/QC#important-parameters>`__),
   you can supply a BED file of regions that you would like to limit the
   operation to. This will do the same thing as in the general
   workaround mentioned above.
-  **computeGCbias**: You can make use of the --filterOut option of
   `computeGCbias <https://github.com/fidelram/deepTools/wiki/All-command-line-options#computegcbias-optional-arguments>`__.
   You will first need to create a BED file that contains all the
   regions you are *not* interested in. Then supply this file of
   RegionsOf\_\_Non\_\_Interest.bed to computeGCbias.

--------------

When should I exclude regions from computeGCbias?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In general, we recommend that you should only correct for GC bias (using
`computeGCbias <https://github.com/fidelram/deepTools/wiki/QC>`__
followed by
`correctGCbias <https://github.com/fidelram/deepTools/wiki/Normalizations>`__)
if you observe that the majority of the genome (the region between
30-60%) is continuously GC-biased **and** you want to compare this
sample with another sample that is not GC-biased.

Sometimes, a certain GC bias is expected, for example for ChIP samples
of H3K4me3 in mammalian samples where GC-rich promoters are expected to
be enriched. To not confound the GC bias caused by the library
preparation with the inherent, expected GC bias, we incorporated the
possibility to supply a file of regions to
`computeGCbias <https://github.com/fidelram/deepTools/wiki/QC>`__ that
will be excluded from the GC bias calculation. This file should
typically contain those regions that one expects to be significantly
enriched per se. This way, the
`computeGCbias <https://github.com/fidelram/deepTools/wiki/QC>`__ will
focus on background regions.

--------------

When should I use bamCoverage, when bamCompare?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Both tools produce bigWig files, i.e. they translate the read-centered
information from a
`BAM <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bam>`__
file into scores for genomic regions of a fixed size. The only
difference is the *number of BAM files* that the tools use as input:
while bamCoverage will only take one BAM file and produce a coverage
file that is mostly normalized for sequencing depth,
`bamCompare <https://github.com/fidelram/deepTools/wiki/Normalizations>`__
will take *two*
`BAM <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bam>`__
files that can be compared with each other using several mathematical
operations.
`bamCompare <https://github.com/fidelram/deepTools/wiki/Normalizations>`__
will always normalize for sequencing depth like bamCoverage, but then it
will perform additional calculations depending on what the user chose,
for example:

-  **bamCompare**:
-  ChIP vs.
   `input <https://github.com/fidelram/deepTools/wiki/Glossary#terminology>`__
   → obtain a bigWig file of log2ratios(ChIP/input)
-  treatment vs. control → obtain a bigWig file of differences
   (Treatment - control)
-  Replicate 1 and Replicate 2 → obtain a bigWig file where the values
   from two BAM files are summed up

--------------

How does computeMatrix handle overlapping genome regions?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| If the
`BED <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bed>`__
file supplied to
`computeMatrix <https://github.com/fidelram/deepTools/wiki/Visualizations>`__
contains regions that overlap, computeMatrix will report those regions
and issue warnings, but they will just be taken as is. If you would like
to avoid that overlapping regions are taken into account, you will need
to clean the
`BED <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bed>`__
file prior to using computeMatrix. There are several possibilities for
modifying your
`BED <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bed>`__
file.
| Let's say your file looks like this:

::

    $ cat testBed.bed
    chr1    10  20  region1
    chr1    7   15  region2
    chr1    18  29  region3
    chr1    35  40  region4
    chr1    10  20  region1Duplicate

-  if you just want to eliminate *identical* entries (here: region1 and
   region1Duplicate), use sort and uniq in the shell (note that the
   label of the identical regions is different - as uniq can only ignore
   fields at the beginning of a file, use rev to revert the sorted file,
   then uniq with ignoring the first field (which is now the name
   column) and then revert back)

   ::

       $ sort -k1,1 -k2,2n testBed.bed | rev | uniq -f1 | rev
       chr1    10  20  region1
       chr1    7   15  region2
       chr1    18  29  region3
       chr1    35  40  region4

-  if you would like to *merge all overlapping regions* to one big one,
   use the BEDtool mergeBed

   -  again, the BED file must be sorted first
   -  -n and -nms tell mergeBed to output the number of overlapping
      regions and the names of them
   -  in the resulting file, regions 1, 2 and 3 are merged

::

    $ sort -k1,1 -k2,2n testBed.bed | mergeBed -i stdin -n -nms 
    chr1    7   29  region2;region1;region1Duplicate;region3    4
    chr1    35  40  region4 1

-  if you would like to *keep only regions that do not overlap* with any
   other region in the same
   `BED <https://github.com/fidelram/deepTools/wiki/Glossary#wiki-bed>`__
   file, use the same mergeBed routine but subsequently filter out those
   regions where several regions were merged

   -  the awk command will check the last field of each line ($NF) and
      will print the original line ($0) only if the last field contained
      a number smaller than 2

::

    $ sort -k1,1 -k2,2n testBed.bed | mergeBed -i stdin -n -nms | awk '$NF < 2 {print $0}'
    chr1    35  40  region4 1

--------------

Why does the maximum value in the heatmap not equal the maximum value in the matrix?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The entire question went like this:

*The max scale in the heatmap does not correspond to the max. values in
the matrix. This suggests that some scaling/normalization is done prior
to heatmap plot. (ex: I have values in the matrix going up to 50 and the
scales on the matrix goes from 0 to ~13)*

Indeed, if you do not indicate --zMax and/or --zMin together with the
heatmapper command, the default of heatmapper is to remove outliers. We
have found that this is beneficial in most cases as outliers tend to
screw up the color schemes. In the code, this looks like this:

::

    if zMin is None:
        matrixFlatten = flattenMatrix(hm.matrixDict)
        # try to avoid outliers by using np.percentile
        zMin = np.percentile(matrixFlatten, 1.0)
        
        if np.isnan(zMin):
          zMin = None

    ...

    if zMax is None:
        if matrixFlatten is None:
          matrixFlatten = flattenMatrix(hm.matrixDict)

    # try to avoid outliers by using np.percentile
       zMax = np.percentile(matrixFlatten, 98.0)

       if np.isnan(zMax):
          zMax = None

If you would like to include the entire range, you can set --zMax and
--zMin to the min and max values you find in your matrix.

--------------

The heatmap I generated looks very "coarse", I would like a much more fine-grained image.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  decrease the
   **`bin <https://github.com/fidelram/deepTools/wiki/Glossary#terminology>`__
   size** when generating the matrix using
   `computeMatrix <https://github.com/fidelram/deepTools/wiki/Visualizations>`__
-  go to "advanced options" → "Length, in base pairs, of the
   non-overlapping
   `bin <https://github.com/fidelram/deepTools/wiki/Glossary#terminology>`__
   for averaging the score over the regions length" → define a smaller
   value, e.g. 50 or 25 bp
-  make sure, however, that you used a sufficiently small
   `bin <https://github.com/fidelram/deepTools/wiki/Glossary#terminology>`__
   size when calculating the bigWig file, though (if generated with
   deepTools, you can check the option
   "`bin <https://github.com/fidelram/deepTools/wiki/Glossary#terminology>`__
   size")

--------------

How can I change the automatic labels of the clusters in a k-means clustered heatmap?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each cluster will get its own box, exactly the same way as different
groups of regions. Therefore, you can use the same option to define the
labels of the final heatmap: In Galaxy: Heatmapper → "Advanced output
options" → "Labels for the regions plotted in the heatmap".

If you indicated 3 clusters for k-means clustering, enter here: C1, C2,
C3 → instead of the full default label ("cluster 1"), the heatmap will
be labeled with the abbreviations.

In the command line, use the --regionsLabel option to define your
customized names.

--------------

How can I manually specify several groups of regions (instead of clustering)?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you would like to compare the profiles and heatmaps for different
groups of regions that you specified yourself, you need to tell
computeMatrix to assess the BED file of regions accordingly. In Galaxy,
you should add one BED file per group (i.e. genes.bed, exons.bed,
introns.bed) while on the command line you should have just one BED file
total with hashed lines indicating the end of a group. For more details,
see the step-by-step-description
`here <https://github.com/fidelram/deepTools/wiki/Visualizations#1st-example-heatmap-with-all-genes-scaled-to-the-one-size-and-user-specified-groups-of-genes>`__.

--------------

What do I have to pay attention to when working with a draft version of a genome?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are working with sequences from a genome that is not included in
our standard descriptions, you need to pay attention to two sets of data
that you might need using deepTools:

#. **Effective genome size** - this is mostly needed for
   `bamCoverage <https://github.com/fidelram/deepTools/wiki/Normalizations>`__
   and
   `bamCompare <https://github.com/fidelram/deepTools/wiki/Normalizations>`__,
   see `below <#effGenomeSize>`__ for details
#. **Reference genome sequence in 2bit format** - this is needed for
   `computeGCbias <https://github.com/fidelram/deepTools/wiki/QC>`__,
   see `below <#2bit>`__ for details

--------------

How do I calculate the effective genome size for an organism that's not in your list?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We plan on including a module that will calculate the effective genome
size for you, but it's not ready yet, so you will have to find a
solution outside of deepTools for the time being.

The "real" effective genome size is the part of the genome that is
*uniquely mappable*. This means that the value will depend on the genome
properties (how many repetitive elements, quality of the assembly etc.)
and the length of the sequenced reads as 100 million 36-bp-reads might
cover less than 100 million 100-bp-reads.

We currently have these options for you:

`1. Use an external tool <#GEM>`__

`2. Use faCount (only if you let reads be aligned non-uniquely,
too!) <#faCount>`__

`3. Use bamCoverage <#mapp_bamCov>`__

`4. Use genomeCoverageBed <#mapp_genomeCov>`__

| 
| **1. Use an external tool**
| There is a tool that promises to calculate the mappability for any
genome given the read length (k-mer length): `**GEM-Mappability
Calculator** <http://algorithms.cnag.cat/wiki/Man:gem-mappability#Mappability.2Falignability>`__.
According to this reply
`here <https://groups.google.com/forum/#!topic/macs-announcement/-iIDkVwenn8>`__,
you can calculate the effective genome size after running this program
by counting the numbers of "!" which stands for uniquely mappable
regions.

| 
| **2. Use faCount**
| If you are using
`bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`__
which reports *multi*-read alignments (= including *non-uniquely*
mapping reads) as a default setting, you can use **faCount from UCSC
tools** to report the total number of bases as well as the number of
bases that are missing from the genome assembly indicated by 'N'. The
effective genome size would then be the total number of base pairs minus
the total number of 'N'.
| Here's an example output of faCount on *D. melanogaster* genome
version dm3:

::

    $ UCSCtools/faCount dm3.fa
    #seq        len     A   C   G    T   N   cpg
    chr2L       23011544    6699731 4811687 4815192  6684734 200     926264
    chr2LHet    368872      90881   58504   57899    90588   71000   10958
    chr2R       21146708    6007371 4576037 4574750  5988450 100     917644
    chr2RHet    3288761     828553  537840   529242  826306  566820  99227
    chr3L       24543557    7113242 5153576  5141498 7135141 100     995078
    chr3LHet    2555491     725986  473888   479000  737434 139183   89647
    chr3R       27905053    7979156 5995211  5980227 7950459 0   1186894
    chr3RHet    2517507     678829  447155   446597  691725  253201  84175
    chr4        1351857     430227  238155   242039  441336  100     43274
    chrU        10049037    2511952 1672330  1672987 2510979 1680789 335241
    chrUextra   29004656    7732998 5109465  5084891 7614402 3462900 986216
    chrX        22422827    6409325 4742952  4748415 6432035 90100   959534
    chrXHet     204112      61961   40017    41813   60321  0    754
    chrYHet     347038      74566   45769    47582   74889  104232   8441
    chrM        19517       8152    2003     1479    7883   0    132
    total       168736537   47352930 33904589 33863611 47246682 6368725 6650479

| In this example:
| Total no. bp = 168,736,537
| Total no. 'N' = 6,368,725

*NOTE*: this method only works if multi-reads are mapped randomly to
their possible locations, in other words if repetitive regions are
covered by reads in which case the effective genome size is the size of
the genome discarding stretches of 'N's.

| 
| **3. Use bamCoverage**
| If you have a sample where you expect the genome to be covered
completely, e.g. from genome sequencing, a very trivial solution is to
use bamCoverage with a bin size of 1 bp and the --outFileFormat option
set to 'bedgraph'. You can then count the number of non-Zero bins (= bp)
which will indicate the mappable genome size for this specific sample.

| 
| **4. Use genomeCoverageBed**
| The BEDtool genomeCoverageBed can be used to calculate the number of
bp in the genome for which 0 reads can be found overlapping. As
described on the `BEDtools
website <http://bedtools.readthedocs.org/en/latest/content/tools/genomecov.html>`__,
you need:

-  a file with the choromosome sizes of your sample's organism
-  a position-sorted BAM file

::

    bedtools genomecov -ibam sortedBAMfile.bam -g genome.size

--------------

Where can I download the 2bit genome files required for *computeGCbias*?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| The 2bit files of most genomes can be found
`here <http://hgdownload.cse.ucsc.edu/gbdb/>`__.
| Search for the .2bit ending. Otherwise, **fasta files can be converted
to 2bit** using the UCSC programm
| faToTwoBit (available for different platforms from
`here <http://hgdownload.cse.ucsc.edu/admin/exe/>`__

--------------

