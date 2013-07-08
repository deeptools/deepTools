from distutils.core import setup
from setuptools import setup

setup(
    name='deepTools',
    version='1.5',
    author='Fidel Ramirez',
    author_email='ramirez@ie-freiburg.mpg.de',
    packages=['deeptools', 'deeptools.test'],
    scripts=['bin/bams2ratio', 'bin/bam2wig', 'bin/correlateBams', 
	'bin/heatmapper', 'bin/bamFingerprint', 'bin/estimateScaleFactor',
	'bin/PE_fragment_size', 'bin/computeMatrix'],
#    url='http://pypi.python.org/pypi/IETools/',
    license='LICENSE.txt',
    description='Useful library to deal with mapped reads in sorted BAM format.',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy >= 1.6.2",
	"pysam >= 0.7.1",
	"bx-python >= 0.5.0",
	"scipy >= 0.10.0"
    ],
)
