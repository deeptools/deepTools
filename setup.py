#-*- coding: utf-8 -*-

import os
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
        print "unable to run git, leaving ecdsa/_version.py alone"
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
        _install.run(self)
#        import ipdb;ipdb.set_trace()
        # check the installation of UCSC tools
        try:
            subprocess.call(["bedGraphToBigWig", "-h"])
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                # handle file not found error.
                import warnings
                # the config file is installed in:
                config_file = self.install_platlib + \
                    "/deeptools/config/deeptools.cfg"
                msg = "\n#############################################\n\n"\
                      "Your computer does not have the UCSC program \n" \
                      "bedGraphToBigWig installed or configured in the \n" \
                      "deepTools config file. In order to output bigwig \n" \
                      "files this program needs to be installed and "\
                      "referred in the \n"\
                      "configuration file located at:\n\n{}\n\n" \
                      "The program can be downloaded from here: " \
                      "http://hgdownload.cse.ucsc.edu/admin/exe/ \n\n" \
                      "The output is set by default to 'bedgraph'\n\n "\
                      "\n#############################################"\
                      "\n\n".format(config_file)
                warnings.warn(msg)
            else:
                # Something else went wrong while
                # trying to run `bedGraphToBigWig`
                raise

setup(
    name='deepTools',
    version=get_version(),
    author='Fidel Ramirez, Friederike Dündar, Björn Grüning, Sarah Diehl',
    author_email='deeptools@googlegroups.com',
    packages=['deeptools', 'deeptools.test'],
    scripts=['bin/bamCompare', 'bin/bamCoverage', 'bin/bamCorrelate',
             'bin/heatmapper', 'bin/bamFingerprint', 'bin/estimateScaleFactor',
             'bin/PE_fragment_size', 'bin/computeMatrix', 'bin/profiler',
             'bin/computeGCBias', 'bin/correctGCBias',
             'bin/bigwigCompare'],
    include_package_data=True,
    package_data={'': ['config/deeptools.cfg']},
#    data_files=[('deepTools', ['./deepTools.cfg'])],
#                ('galaxy', ['galaxy/bamCompare.xml','galaxy/bamCoverage.xml',
#                            'galaxy/heatmapper.xml'])],
    url='http://pypi.python.org/pypi/deepTools/',
    license='LICENSE.txt',
    description='Useful library to deal with mapped reads in sorted '
    'BAM format.',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy >= 1.6.2",
        "scipy >= 0.10.0",
        "matplotlib >= 1.2",
        "pysam >= 0.7.1",
        "bx-python >= 0.5.0",
    ],
    cmdclass={'sdist': sdist, 'install': install}
)
