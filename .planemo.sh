#!/bin/bash
blah=`mktemp -d`
/home/travis/build/fidelram/deepTools/foo/bin/planemo conda_init --conda_prefix $blah/conda
export PATH=$blah/conda/bin:$PATH
conda create -y --name deeptools_galaxy numpy matplotlib scipy
source activate deeptools_galaxy
conda install -c bioconda samtools
git clone --depth 1 https://github.com/galaxyproject/galaxy.git clone
cd clone
./scripts/common_startup.sh --skip-venv --dev-wheels
cd ..
pip install . 
/home/travis/build/fidelram/deepTools/foo/bin/planemo test --galaxy_root clone --test_data galaxy/wrapper/test-data/ --skip_venv \
galaxy/wrapper/bamCompare.xml \
galaxy/wrapper/bamCoverage.xml \
galaxy/wrapper/bamPEFragmentSize.xml \
galaxy/wrapper/bigwigCompare.xml \
galaxy/wrapper/computeGCBias.xml \
galaxy/wrapper/computeMatrix.xml \
galaxy/wrapper/correctGCBias.xml \
galaxy/wrapper/multiBamSummary.xml \
galaxy/wrapper/multiBigwigSummary.xml \
galaxy/wrapper/plotCoverage.xml \
galaxy/wrapper/plotFingerprint.xml 2>&1 | grep -v -e "^galaxy" | grep -v -e "^requests"
