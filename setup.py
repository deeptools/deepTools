# -*- coding: utf-8 -*-

import os
import subprocess
import re
from distutils import sysconfig
import glob

from setuptools import setup, Extension, find_packages
from setuptools.command.sdist import sdist as _sdist
from setuptools.command.install import install as _install

VERSION_PY = """
# This file is originally generated from Git information by running 'setup.py
# version'. Distribution tarballs contain a pre-generated copy of this file.

__version__ = '%s'
"""

srcs = [x for x in glob.glob("deeptoolsintervals/tree/*.c")]

libs = ["z"]
additional_libs = [sysconfig.get_config_var("LIBDIR"), sysconfig.get_config_var("LIBPL")]

module1 = Extension('deeptoolsintervals.tree',
                    sources=srcs,
                    libraries=libs,
                    library_dirs=additional_libs,
                    include_dirs=[sysconfig.get_config_var("INCLUDEPY")])


def update_version_py():
    if not os.path.isdir(".git"):
        print("This does not appear to be a Git repository.")
        return
    try:
        p = subprocess.Popen(["git", "describe",
                              "--tags", "--always"],
                             stdout=subprocess.PIPE)
    except EnvironmentError:
        print("unable to run git, leaving deeptools/_version.py alone")
        return
    stdout = p.communicate()[0]
    if p.returncode != 0:
        print("unable to run git, leaving deeptools/_version.py alone")
        return
    ver = stdout.decode().strip().replace("-g", "-")
    f = open("deeptools/_version.py", "w")
    f.write(VERSION_PY % ver)
    f.close()
    print("set deeptools/_version.py to '%s'" % ver)


def get_version():
    try:
        f = open("deeptools/_version.py")
    except EnvironmentError:
        return None
    for line in f.readlines():
        mo = re.match("__version__ = '([^']+)'", line)
        if mo:
            ver = mo.group(1)
            return ver
    return None


class sdist(_sdist):

    def run(self):
        update_version_py()
        self.distribution.metadata.version = get_version()
        return _sdist.run(self)


class install(_install):

    def run(self):
        update_version_py()
        self.distribution.metadata.version = get_version()
        _install.run(self)
        return


def openREADME():
    """
    This is only needed because README.rst is UTF-8 encoded and that won't work
    under python3 iff sys.getfilesystemencoding() returns 'ascii'

    Since open() doesn't accept an encoding in python2...
    """
    try:
        f = open("README.rst", encoding="utf-8")
    except:
        f = open("README.rst")

    foo = f.read()
    f.close()
    return foo


setup(
    name='deepTools',
    version=get_version(),
    author='Fidel Ramirez,  Devon P Ryan, Björn Grüning, Friederike Dündar, Sarah Diehl,'
    ' Vivek Bhardwaj, Fabian Kilpert, Andreas S Richter, Steffen Heyne, Thomas Manke',
    author_email='deeptools@googlegroups.com',
    packages=find_packages(),
    scripts=['bin/bamCompare', 'bin/bamCoverage', 'bin/multiBamSummary',
             'bin/plotHeatmap', 'bin/plotFingerprint', 'bin/estimateScaleFactor',
             'bin/bamPEFragmentSize', 'bin/computeMatrix', 'bin/plotProfile',
             'bin/computeGCBias', 'bin/correctGCBias', 'bin/multiBigwigSummary',
             'bin/bigwigCompare', 'bin/plotCoverage', 'bin/plotPCA', 'bin/plotCorrelation',
             'bin/plotEnrichment', 'bin/deeptools', 'bin/computeMatrixOperations',
             'bin/estimateReadFiltering', 'bin/alignmentSieve'],
    include_package_data=True,
    url='http://pypi.python.org/pypi/deepTools/',
    license='LICENSE.txt',
    description='Useful tools for exploring deep sequencing data ',
    long_description=openREADME(),
    classifiers=[
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    install_requires=[
        "numpy >= 1.9.0",
        "scipy >= 0.17.0",
        "matplotlib >= 2.1.2",
        "pysam >= 0.14.0",
        "numpydoc >=0.5",
        "pyBigWig >=0.2.1",
        "py2bit >= 0.2.0",
        "plotly >= 1.9.0"
    ],
    zip_safe=False,
    ext_modules=[module1],
    cmdclass={'sdist': sdist, 'install': install}
)
