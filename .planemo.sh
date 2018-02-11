#!/bin/bash
echo $PATH
planemo=./miniconda/envs/planemo/bin/planemo
owd=`pwd`
temp_dir=`mktemp -d`
cd $temp_dir
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -y --name gxtest numpy pysam
#source activate gxtest
#git clone --depth 1 https://github.com/galaxyproject/galaxy.git
#cd galaxy
#make client
#./scripts/common_startup.sh --skip-venv --dev-wheels
#cd ..
## reset what's available in conda
#cd $owd
#conda install --yes -c conda-forge numpy scipy matplotlib==2.1.0 plotly==2.0.12 cython
#conda install --yes -c bioconda -c conda-forge pysam==0.14.0 pyBigWig py2bit
#python setup.py install

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
