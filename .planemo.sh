#!/bin/bash
# Some versions of planemo don't handle symlinks
unlink galaxy/wrapper/test-data/test.bw
cp deeptools/test/test_heatmapper/test.bw galaxy/wrapper/test-data/test.bw

if [[ $1 == "1" ]] ; then
    wrappers="galaxy/wrapper/alignmentSieve.xml \
    galaxy/wrapper/bamCompare.xml \
    galaxy/wrapper/bamCoverage.xml \
    galaxy/wrapper/bamPEFragmentSize.xml \
    galaxy/wrapper/bigwigCompare.xml \
    galaxy/wrapper/computeGCBias.xml"
elif [[ $1 == "2" ]] ; then
    wrappers="galaxy/wrapper/computeMatrix.xml \
    galaxy/wrapper/computeMatrixOperations.xml \
    galaxy/wrapper/correctGCBias.xml \
    galaxy/wrapper/estimateReadFiltering.xml \
    galaxy/wrapper/multiBamSummary.xml \
    galaxy/wrapper/multiBigwigSummary.xml"
else
    wrappers="galaxy/wrapper/plotCorrelation.xml \
    galaxy/wrapper/plotCoverage.xml \
    galaxy/wrapper/plotEnrichment.xml \
    galaxy/wrapper/plotFingerprint.xml \
    galaxy/wrapper/plotHeatmap.xml \
    galaxy/wrapper/plotPCA.xml \
    galaxy/wrapper/plotProfiler.xml"
fi

planemo lint ${wrappers}
planemo test --no_dependency_resolution --galaxy_branch release_20.09 --install_galaxy ${wrappers} 2>&1 | grep -v -e "^galaxy" | grep -v -e "^requests"
test ${PIPESTATUS[0]} -eq 0
