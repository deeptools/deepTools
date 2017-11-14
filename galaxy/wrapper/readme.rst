========================
Galaxy deeptools wrapper
========================

deepTools are user-friendly tools for the normalization and visualization of 
deep-sequencing data.
They address the challenge of visualizing the large amounts of data that are now
routinely generated from sequencing centers in a meaningful way. 
To do so, deepTools contain useful routines to process the mapped reads data 
through removal of duplicates and different filtering options to create coverage
files in standard bedGraph and bigWig file formats. deepTools allow the creation
of normalized coverage files or the comparison between two files 
(for example, treatment and control). Finally, using such normalized and 
standardized files, multiple visualizations can be created to identify 
enrichments with functional annotations of the genome. 
For a gallery of images that can be produced and a description 
of the tools see our poster_.

.. _poster: http://f1000.com/posters/browse/summary/1094053

deeptools is developed under here:

    https://github.com/deeptools/deepTools

For support, questions, or feature requests contact: deeptools@googlegroups.com


============
Installation
============

Requirements: python-2.7

Galaxy should be able to automatically install all other dependencies, such as numpy or scipy.

For the best performance we recommend to install blas/lapack/atlas in your environment before
installing deepTools from the Tool Shed.


========
Citation
========

deeptools are currently under review. In the meantime please refere to https://github.com/deeptools/deepTools.


=======
History
=======

 * v1.0:        Initial public release
 * v1.5.8.2:    Include new citation tag, update version to 1.5.8.2 and change wrapper version


Licence (MIT)
=============

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
