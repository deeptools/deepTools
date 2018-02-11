#!/bin/bash
planemo=$HOME/miniconda/envs/planemo/bin/planemo

#galaxy/wrapper/correctGCBias.xml \
$planemo test --no_dependency_resolution --install_galaxy \
galaxy/wrapper/alignmentSieve.xml \
galaxy/wrapper/bamCompare.xml \
galaxy/wrapper/bamCoverage.xml \
galaxy/wrapper/bamPEFragmentSize.xml \
galaxy/wrapper/bigwigCompare.xml \
galaxy/wrapper/computeGCBias.xml \
galaxy/wrapper/computeMatrix.xml \
galaxy/wrapper/computeMatrixOperations.xml \
galaxy/wrapper/estimateReadFiltering.xml \
galaxy/wrapper/multiBamSummary.xml \
galaxy/wrapper/multiBigwigSummary.xml \
galaxy/wrapper/plotCorrelation.xml \
galaxy/wrapper/plotCoverage.xml \
galaxy/wrapper/plotEnrichment.xml \
galaxy/wrapper/plotFingerprint.xml \
galaxy/wrapper/plotHeatmap.xml \
galaxy/wrapper/plotPCA.xml \
galaxy/wrapper/plotProfiler.xml 2>&1 | grep -v -e "^galaxy" | grep -v -e "^requests"
test ${PIPESTATUS[0]} -eq 0
