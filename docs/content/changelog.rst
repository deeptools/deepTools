Changes in deepTools4.0
=======================

# Changes

## bamCoverage

 - --no_collapse flag to not merge bins with equal coverage values together.
 
## computeMatrix

 - --sortRegions 'no' option no longer exists
 - Sorting ascend / descend no longer has subsorting by position.
 - --quiet / -q option no longer exists.
 - bed files in computeMatrix no longer support '#' to define groups.
 - 'chromosome matching' i.e. chr1 <-> 1, chrMT <-> MT is no longer performed.

## normalization

 - Exactscaling is no longer an option, it's always performed.

## alignmentSieve

- options label, smartLabels, genomeChunkLength are removed.
- ignoreDuplicates is removed, and (if wanted) should be set by the SamFlagExclude setting.












# Good to know
















# Testing
 
## computeMatrix
 - referencePoint: TSS, center, TES
 - sortRegions: descend, ascend, keep
 - sortUsing: mean, median, max, min, sum, region_length
 - averageTypeBins: mean, median, min, max ,std, sum
 - skipZeros
 - duplicate renaming _r1, _r2, ...
 - GTF, BED3, BED6, BED12, mixedBED (?)
 - scaleRegions, un5, un3, regionbodylength, metagene

 ## alignmentSieve

 - unmapped reads to unfiltered_out

## multiBamSummary
 - npz output has labels encoded as u8s, no longer strings.

# Todo

- AlignmentSieve: Shift, Bed, Optimization.
- bamCoverage / bamCompare: filtering, extend.