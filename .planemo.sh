#!/bin/bash
blah=`mktemp -d`
conda config --add channels conda-forge
conda config --add channels bioconda
conda install -c bioconda samtools
/home/travis/build/fidelram/deepTools/foo/bin/planemo database_create galaxy

git clone --depth 1 --single-branch --branch release_16.10 https://github.com/galaxyproject/galaxy.git clone
cd clone
touch tool-data/twobit.loc
#Add the custom data types
sed -i '4i\    <datatype extension="deeptools_compute_matrix_archive" type="galaxy.datatypes.binary:CompressedArchive" subclass="True" display_in_upload="True"/>' config/datatypes_conf.xml.sample
sed -i '5i\    <datatype extension="deeptools_coverage_matrix" type="galaxy.datatypes.binary:CompressedArchive" subclass="True" display_in_upload="True"/>' config/datatypes_conf.xml.sample
./scripts/common_startup.sh --skip-venv --dev-wheels
cd ..
pip install .
/home/travis/build/fidelram/deepTools/foo/bin/planemo test --galaxy_root clone --test_data galaxy/wrapper/test-data/ --skip_venv --postgres --no_conda_auto_install --no_conda_auto_init \
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
