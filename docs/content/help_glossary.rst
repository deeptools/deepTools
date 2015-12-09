Glossary
========

Like most specialized fields, next-generation sequencing has inspired many an acronym.
We are trying to keep track of those :ref:`abbreviations` that we heavily use.
Do make us aware if something is unclear: deeptools@googlegroups.com

We have also assembled a short glossary of common :ref:`terminology` and
:ref:`file formats` of next-generation sequencing data.

.. _abbreviations:

Abbreviations
---------------

Reference genomes are usually referred to by their abbreviations, such as:

* hg19 = human genome, version 19
* mm9 = *Mus musculus* genome, version 9
* dm3 = *Drosophila melanogaster*, version 3
* ce10 = *Caenorhabditis elegans*, version 10

For a more comprehensive list of available reference genomes and their abbreviations,
see the `UCSC data base <http://hgdownload.soe.ucsc.edu/downloads.html>`_.
 
+---------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| Acronym       | full phrase                              | Synonyms/Explanation                                                                                                                          |
+===============+==========================================+===============================================================================================================================================+
| <ANYTHING>-seq| -sequencing                              |indicates that an experiment was completed by DNA sequencing using NGS                                                                         |
+---------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| ChIP-seq      | chromatin immunoprecipitation sequencing | NGS technique for detecting transcription factor binding sites and histone modifications (see entry *Input* for more information)             |
+---------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| DNase         | deoxyribonuclease I                      | DNase I digestion is used to determine active ("open") chromatin regions                                                                      |
+---------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| HTS           | high-throughput sequencing               | next-generation sequencing, massive parallel short read sequencing, deep sequencing                                                           |
+---------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| MNase         | micrococcal nuclease                     | MNase digestion is used to determine sites with nucleosomes                                                                                   |
+---------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| NGS           | next-generation sequencing               | high-throughput (DNA) sequencing, massive parallel short read sequencing, deep sequencing                                                     |
+---------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| RPGC          | reads per genomic content                | normalize reads to 1x sequencing depth, sequencing depth is defined as: (mapped reads x fragment length) / effective genome size              |
+---------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+
| RPKM          | reads per kilobase per million reads     | normalize read numbers: RPKM (per bin) = reads per bin / ( mapped reads (in millions) x bin length (kb))                                      |
+---------------+------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------+

For a review of popular *-seq applications, see `Zentner and Henikoff`_.

.. _terminology:

NGS terminology
---------------
In addition to abbreviations, there are many specialized terms and different labs and people use different terms for the same thing.

**bin**

* synonyms: window, region
* For many calculations, the genome is divided into smaller chunks, for example for the calculation of read coverages. These regions can be as small as 1 bp, but really, they can be any size. We commonly use the term "bin" for those artificially created genome parts as we feel that we "store" scores (e.g. read coverages or motif scores) in them.


**Input**

* Control experiment typically done for ChIP-seq experiments 
* While ChIP-seq relies on antibodies to enrich for DNA fragments bound to a certain protein, the input sample should be processed exactly the same way, excluding the antibody. This way, one hopes to account for biases introduced by the sample handling and the general chromatin structure of the cells 

**read**

* synonym: tag
* This term refers to the piece of DNA that is sequenced ("read") by the sequencers. We try to differentiate between "read" and "DNA fragment" as the fragments that are put into the sequencer tend to be in the range of 200-1000 bp of which only the first 30 to 100 bp (depending on the read length) are in fact sequenced. Most of the deepTools will not only take those 30 to 100 bp into consideration when calculating coverages, instead they will extend the reads to match the original DNA fragment size. (The original size will either be given by you or, if you used paired-end sequencing, can be calculated by the distance of two read mates).

.. _file formats:

File Formats
-------------------

Data obtained from next-generation sequencing data must be processed several times.
Most of the processing steps are aimed at extracting only those information that are
truly needed for a specific down-stream analysis and to discard all the redundant entries.
Therefore, **specific data formats are often associated with different steps of a data processing pipeline**.
These associations, are by no means binding, but you should understand what kind of data is represented in which data format
- this will help you to select the correct tools further down the road.

Here, we just want to give very brief key descriptions of the file, for elaborate information we will link to external websites.
Be aware, that the file name sorting here is purely alphabetically, not according to their usage within an analysis pipeline that is depicted here:

.. image:: images/flowChart_FileFormats.png

Follow the links for more information on the different tool collections mentioned in the figure:

`samtools <http://samtools.sourceforge.net/http://samtools.sourceforge.net/>`_ |
`UCSCtools <http://hgdownload.cse.ucsc.edu/admin/exe/>`_ |
`BEDtools <http://bedtools.readthedocs.org/en/latest/>`_ |

.. _2bit:

2bit
^^^^

* compressed, binary version of genome sequences that are often stored in :ref:`FASTA` 
* most genomes in 2bit format can be found `at UCSC <http://hgdownload.cse.ucsc.edu/gbdb/>`_
* :ref:`FASTA` files can be converted to 2bit using the UCSC programm *faToTwoBit*, which is available for different platforms at `UCSC <http://hgdownload.cse.ucsc.edu/admin/exe/>`_
* more information can be found `here <http://jcomeau.freeshell.org/www/genome/2bitformat.html>`_ or from `UCSC <http://genome.ucsc.edu/FAQ/FAQformat.html#format7>`_

.. _BAM:

BAM
^^^^

* typical file extension: .bam
* *binary* file format (complement to :ref:`SAM`)
* contains information about sequenced reads *after alignment* to a reference genome
* each line = 1 mapped read, with information about:
    *  its mapping quality (how certain is the read alignment to this particular genome locus?)
    *  its sequencing quality (how well was each base pair detected during sequencing?)
    *  its DNA sequence
    *  its location in the genome
    *  etc.
* highly recommended format for storing data
* to make a BAM file human-readable, one can, for example, use the program *samtools view* 
* for more information, see below for the definition of :ref:`SAM` files

.. _bed:

bed
^^^

* typical file extension: .bed
* text file
* used for genomic intervals, e.g. genes, peak regions etc.
* actually, there is a rather strict definition of the format that can be found at `UCSC`_
* for deepTools, the first 3 columns are important: chromosome, start position of the region, end position of the genome
* do not confuse it with the :ref:`bedgraph` format (eventhough they are quite similar)
* example lines from a BED file of mouse genes (note that the start position is 0-based, the end-position 1-based, following UCSC conventions for BED files):
::

    chr1    3204562 3661579 NM_001011874 Xkr4   -
    chr1    4481008 4486494 NM_011441    Sox17  -
    chr1    4763278 4775807 NM_001177658 Mrpl15 -
    chr1    4797973 4836816 NM_008866    Lypla1 +


.. _bedgraph:

bedGraph 
^^^^^^^^

* typical file extension: .bg, .bedgraph
* text file
* similar to BED file (not the same!), it can *only* contain 4 columns and the 4th column *must* be a score
* again, read the `UCSC description <https://genome.ucsc.edu/FAQ/FAQformat.html#format1.8>`_  for more details
* 4  exemplary lines from a bedGraph file (like BED files following the UCSC convention, the start position is 0-based, the end-position 1-based in bedGraph files):
::

    chr1 10 20 1.5
    chr1 20 30 1.7
    chr1 30 40 2.0
    chr1 40 50 1.8

.. _bigwig:

bigWig 
^^^^^^

* typical file extension: .bw, .bigwig
* *binary* version of a :ref:`bedgraph` file
* usually contains 4 columns: chromosome, start of genomic bin, end of genomic bin, score
* the score can be anything, e.g. an average read coverage
* `UCSC description <https://genome.ucsc.edu/FAQ/FAQformat.html#format6.1>`_ for more details

.. _FASTA:

FASTA 
^^^^^^

* typical file extension: .fasta
* text file, often gzipped (.fasta.gz)
* very simple format for **DNA/RNA** or **protein** sequences, this can be anything from small pieces of DNA or proteins to entire genome information (most likely, you will get the genome sequence of your organism of interest in fasta format)
* see the :ref:`2bit` file format entry for a compressed alternative of the fasta format
* example from [wikipedia](http://en.wikipedia.org/wiki/FASTA_format "wikipedia entry on FASTA files") showing exactly one sequence:
::

    >gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
     LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV
     EWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG
     LLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVIL
     GLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX
     IENY

.. _fastq:

FASTQ
^^^^^

* typical file extension: .fastq, fq
* text file, often gzipped (--> .fastq.gz)
* contains raw read information -- 4 lines per read:
	 * read ID
	 * base calls
	 * additional information or empty line
	 * sequencing quality measures - 1 per base call
* note that there are no information about where in the genome the read originated from
* example from the `wikipedia page <http://en.wikipedia.org/wiki/Fastq>`_
::

  # a FASTQ file containing a single sequence might look like this:
    
    @read001														# read ID
    GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT	# read sequence
    +																# usually empty line
    !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65	# ASCII-encoded quality scores

* the quality string (here, the line starting with !) contains `ASCII-encoded <http://www.asciitable.com/>`_ measures of probability that the respective base call was correct
	* the ASCII-encoding usually starts with an offset of at least 33 because the first 32 ASCII representations correspond to invisible entries, such as white spaces
	* thus, a base call quality of 0 will be 0 + 33 = 33; in ASCII, 33 corresponds to the ! sign -- the first base in the above shown example is therefore a base where the base identification was rather uncertain
	* however, over the years, Illumina has changed the offset a couple of times, which may sometimes lead to problems with aligners
	
.. image:: images/glossary_ascii.png
	
* if you need to find out what type of ASCII-encoding your .fastq file contains, you can simply run `FastQC <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ -- its summery file will tell you

.. _SAM:

SAM
^^^ 

* typical file extension: .sam
* usually the result of an alignment of sequenced reads to a reference genome
* contains a short header section (entries are marked by @ signs) and an alignment section where each line corresponds to a single read (thus, there can be millions of these lines)

.. image:: images/glossary_sam.png

* **header section**:
	* tab-delimited lines, beginning with @, followed by tag:value pairs
	* *tag* = two-letter string that defines the content and the format of *value*
	
* **alignment section**:
	* each line contains information about its mapping quality, its sequence, its location in the genome etc.
::

    r001 163 chr1 7 30 8M2I4M1D3M = 37 39 TTAGATAAAGGATACTG *
    r002 0 chr1 9 30 3S6M1P1I4M * 0 0 AAAAGATAAGGATA *

	* the **flag in the second field** contains the answer to several yes/no assessments that are encoded in a single number
	* for more details on the flag, see `this thorough explanation <http://ppotato.wordpress.com/2010/08/25/samtool-bitwise-flag-paired-reads/>`_ or `this more technical explanation <http://blog.nextgenetics.net/?e=18>`_
    * the **CIGAR string in the 6th field** represents the types of operations that were needed in order to align the read to the specific genome location:
		* insertion
		* deletion (small deletions denoted with `D`, bigger deletions, e.g., for spliced reads, denoted with `N`)
		* clipping (deletion at the ends of a read)

.. warning::  Although the SAM/BAM format is rather meticulously defined and documented, whether an alignment program 
will produce a SAM/BAM file that adheres to these principles is completely up to the programmer. The mapping score, CIGAR string,
and particularly, all optional flags (fields >11) are often very differently defined depending on the program. If you plan on filtering
your data based on any of these criteria, make sure you know exactly how these entries were calculated!

.. _UCSC: <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>
.. _Zentner and Henikoff: <http://genomebiology.com/2012/13/10/250>