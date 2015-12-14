# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import re

from setuptools import setup
from setuptools.command.sdist import sdist as _sdist
from setuptools.command.install import install as _install

VERSION_PY = """
# This file is originally generated from Git information by running 'setup.py
# version'. Distribution tarballs contain a pre-generated copy of this file.

__version__ = '%s'
"""


def update_version_py():
    if not os.path.isdir(".git"):
        print "This does not appear to be a Git repository."
        return
    try:
        p = subprocess.Popen(["git", "describe",
                              "--tags", "--always"],
                             stdout=subprocess.PIPE)
    except EnvironmentError:
        print "unable to run git, leaving deeptools/_version.py alone"
        return
    stdout = p.communicate()[0]
    if p.returncode != 0:
        print "unable to run git, leaving deeptools/_version.py alone"
        return
    ver = stdout.strip()
    f = open("deeptools/_version.py", "w")
    f.write(VERSION_PY % ver)
    f.close()
    print "set deeptools/_version.py to '%s'" % ver


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

        # check installation of several components
        samtools_installed = self.checkProgramIsInstalled(
            'samtools', 'view',
            'http://samtools.sourceforge.net/',
            'correctGCbias')

        bedGraphToBigWig_installed = self.checkProgramIsInstalled(
            'bedGraphToBigWig', '-h',
            'http://hgdownload.cse.ucsc.edu/admin/exe/',
            'bamCoverage, bamCompare, correctGCbias')

        if not samtools_installed or not bedGraphToBigWig_installed:
            msg = "\n##########\nSome tools were not fund.\n"\
                "If you already have a copy of these programs installed\n"\
                "please be sure that they are found in your PATH or\n"\
                "that they are referred to in the configuration file of deepTools\n"\
                "located at:\n\n {}\n\n".format(self.config_file)
            sys.stderr.write(msg)

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

setup(
    name='deepTools',
    version=get_version(),
    author='Fidel Ramirez, Friederike Dündar, Björn Grüning, Sarah Diehl',
    author_email='deeptools@googlegroups.com',
    packages=['deeptools', 'deeptools/config', 'deeptools.test'],
    scripts=['bin/bamCompare', 'bin/bamCoverage', 'bin/bamCorrelate',
             'bin/plotHeatmap', 'bin/bamFingerprint', 'bin/estimateScaleFactor',
             'bin/bamPEFragmentSize', 'bin/computeMatrix', 'bin/plotProfile',
             'bin/computeGCBias', 'bin/correctGCBias', 'bin/bigwigCorrelate',
             'bin/bigwigCompare', 'bin/plotCoverage', 'bin/plotPCA', 'bin/plotCorrelation'],
    include_package_data=True,
    package_data={'': ['config/deeptools.cfg']},
    url='http://pypi.python.org/pypi/deepTools/',
    license='LICENSE.txt',
    description='Useful library to deal with mapped reads in sorted '
    'BAM format.',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy >= 1.8.0",
        "scipy >= 0.15.0",
        "matplotlib >= 1.4.0",
        "pysam >= 0.8.2",
        "bx-python >= 0.5.0",
        "numpydoc >=0.5",
        "pyBigWig >=0.1.9"
    ],
    cmdclass={'sdist': sdist, 'install': install}
)
