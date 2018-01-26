#!/bin/bash
#!/bin/bash
planemo=./foo/bin/planemo
owd=`pwd`
temp_dir=`mktemp -d`
cd $temp_dir
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -y --name gxtest numpy bx-python pysam
source activate gxtest
git clone --depth 1 https://github.com/galaxyproject/galaxy.git
cd galaxy
make client
./scripts/common_startup.sh --skip-venv --dev-wheels
cd ..
# reset what's available in conda
cd $owd
conda install --yes -c conda-forge python=2.7 numpy scipy matplotlib==2.1.0 nose flake8 plotly==2.0.12
conda install --yes -c bioconda -c conda-forge pysam pyBigWig py2bit planemo
python setup.py install

#galaxy/wrapper/correctGCBias.xml \
/home/travis/build/deeptools/deepTools/foo/bin/planemo test --galaxy_root clone --test_data galaxy/wrapper/test-data/ --skip_venv --postgres --no_conda_auto_install --no_conda_auto_init \
galaxy/wrapper/bamCompare.xml \
galaxy/wrapper/bamCoverage.xml \
galaxy/wrapper/bamPEFragmentSize.xml \
galaxy/wrapper/bigwigCompare.xml \
galaxy/wrapper/computeGCBias.xml \
galaxy/wrapper/computeMatrix.xml \
galaxy/wrapper/computeMatrixOperations.xml \
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
