#!/bin/bash
blah=`mktemp -d`
/home/travis/build/fidelram/deepTools/foo/bin/planemo database_create galaxy
bash miniconda.sh -b -p $blah/conda
export PATH=$blah/conda/bin:$PATH
conda create -y --name deeptools_galaxy numpy matplotlib scipy
source activate deeptools_galaxy
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda samtools
git clone --depth 1 --single-branch --branch release_16.10 https://github.com/galaxyproject/galaxy.git clone
cd clone
#Add the custom data types
sed -i '4i\    <datatype extension="deeptools_compute_matrix_archive" type="galaxy.datatypes.binary:CompressedArchive" subclass="True" display_in_upload="True"/>' config/datatypes_conf.xml.sample
sed -i '5i\    <datatype extension="deeptools_coverage_matrix" type="galaxy.datatypes.binary:CompressedArchive" subclass="True" display_in_upload="True"/>' config/datatypes_conf.xml.sample
pip install PyYAML==3.11
./scripts/common_startup.sh --skip-venv --dev-wheels
cd ..
pip install .
/home/travis/build/fidelram/deepTools/foo/bin/planemo test --galaxy_root clone --test_data galaxy/wrapper/test-data/ --skip_venv --postgres --conda_prefix $blah/conda --conda_ensure_channels bioconda,conda-forge \
galaxy/wrapper/bamCompare.xml \
galaxy/wrapper/bamCoverage.xml \
galaxy/wrapper/bamPEFragmentSize.xml \
galaxy/wrapper/bigwigCompare.xml \
galaxy/wrapper/computeGCBias.xml \
galaxy/wrapper/computeMatrix.xml \
galaxy/wrapper/computeMatrixOperations.xml \
galaxy/wrapper/correctGCBias.xml \
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
