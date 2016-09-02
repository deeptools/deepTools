# -*- coding: utf-8 -*-

import os
import sys
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
        if os.environ.get('DEEP_TOOLS_NO_CONFIG', False):
            return
        self.config_file = self.install_platlib + \
            "/deeptools/config/deeptools.cfg"

    def checkProgramIsInstalled(self, program, args, where_to_download,
                                affected_tools):
        try:
            subprocess.Popen([program, args],
                             stderr=subprocess.PIPE,
                             stdout=subprocess.PIPE)
            return True
        except EnvironmentError:
            # handle file not found error.
            # the config file is installed in:
            msg = "\n**{0} not found. This " \
                  "program is needed for the following "\
                  "tools to work properly:\n"\
                  " {1}\n"\
                  "{0} can be downloaded from here:\n " \
                  " {2}\n".format(program, affected_tools,
                                  where_to_download)
            sys.stderr.write(msg)

        except Exception as e:
            sys.stderr.write("Error: {}".format(e))


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
             'bin/plotEnrichment', 'bin/deeptools'],
    include_package_data=True,
    package_data={'': ['config/deeptools.cfg']},
    url='http://pypi.python.org/pypi/deepTools/',
    license='LICENSE.txt',
    description='Useful tools for exploring deep sequencing data ',
    long_description=openREADME(),
    install_requires=[
        "numpy >= 1.9.0",
        "scipy >= 0.17.0",
        "matplotlib >= 1.4.0",
        "pysam >= 0.8.2",
        "twobitreader >= 3.1.4",
        "numpydoc >=0.5",
        "pyBigWig >=0.2.1"
    ],
    zip_safe=False,
    ext_modules=[module1],
    cmdclass={'sdist': sdist, 'install': install}
)
