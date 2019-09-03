Effective Genome Size
=====================

A number of tools can accept an "effective genome size". This is defined as the length of the "mappable" genome. There are two common alternative ways to calculate this::

  1. The number of non-N bases in the genome.
  2. The number of regions (of some size) in the genome that are uniquely mappable (possibly given some maximal edit distance).

Option 1 can be computed using ``faCount`` from `Kent's tools <http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/>`__. The effective genome size for a number of genomes using this method is given below:

======== ==============
Genome   Effective size
======== ==============
GRCh37   2864785220
GRCh38   2913022398
GRCm37   2620345972
GRCm38   2652783500
dm3      162367812
dm6      142573017
GRCz10   1369631918
WBcel235 100286401
TAIR10   119481543 
======== ==============

These values only appropriate if multimapping reads are included. If they are excluded (or there's any MAPQ filter applied), then values derived from option 2 are more appropriate. These are then based on the read length. We can approximate these values for various read lengths using the `khmer program <http://khmer.readthedocs.io/en/v2.1.1/>`__ program and ``unique-kmers.py`` in particular. A table of effective genome sizes given a read length using this method is provided below:

=========== ========== ========== ========== ========== ========= ========= ========== ========
Read length GRCh37     GRCh38     GRCm37     GRCm38     dm3       dm6       GRCz10     WBcel235
=========== ========== ========== ========== ========== ========= ========= ========== ========
50          2685511504 2701495761 2304947926 2308125349 130428560 125464728 1195445591 95159452
75          2736124973 2747877777 2404646224 2407883318 135004462 127324632 1251132686 96945445
100         2776919808 2805636331 2462481010 2467481108 139647232 129789873 1280189044 98259998
150         2827437033 2862010578 2489384235 2494787188 144307808 129941135 1312207169 98721253
200         2855464000 2887553303 2513019276 2520869189 148524010 132509163 1321355241 98672758
=========== ========== ========== ========== ========== ========= ========= ========== ========
