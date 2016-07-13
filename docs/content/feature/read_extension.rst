Read extension
==============

In the majority of NGS experiment, DNA (or RNA) is fragmented into small stretches and only the ends of these fragments sequenced. For many applications, it's desirable to quantify coverage of the entire original fragments over the genome. Consequently, there is an `--extendReads` option present throughout deepTools. This works as follows:

Paired-end reads
----------------

 1. Regions of the genome are sampled to determine the median fragment/read length.
 2. The genome is subdivided into disjoint regions. Each of these regions comprises one or more bins of some desired size (specified by `-bs`).
 3. For each region, all alignments overlapping it are gathered. In addition, all alignments within 2000 bases are gathered, as 2000 bases is the maximum allowed fragment size.
 4. The resulting collection of alignments are all extended according to their fragment length, which for paired-end reads is indicated in BAM files.

   - For singletons, the expected fragment length from step 1 is used.

 5. For each of the extended reads, the count in each bin that it overlaps is incremented.

Single-end reads
----------------

 1. An extension length, L, is specified.
 2. The genome is subdivided into disjoint regions. Each of these regions comprises one or more bins of some desired size (specified by `-bs`).
 3. For each region, all alignments overlapping it are gathered. In addition, all alignments within 2000 bases are gathered, as 2000 bases is the maximum allowed fragment size.
 4. The resulting collection of alignments are all extended to length L.
 5. For each of the extended reads, the count in each bin that it overlaps is incremented.

Blacklisted regions
-------------------

The question likely arises as to how alignments originating inside of blacklisted regions are handled. In short, any alignment contained completely within a blacklisted region is ignored, regardless of whether it would extend into a non-blacklisted region or not. Alignments only partially overlapping blacklisted regions are treated as normal, as are pairs of reads that span over a blacklisted region. This is primarily for the sake of performance, as otherwise each extended read would need to be checked to see if it overlaps a blacklisted region.
