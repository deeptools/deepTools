Accessing datasets hosted on deepBlue
=====================================

`deepBlue <http://dx.doi.org/10.1093/nar/gkw211>`__ is an epigenome dataset server hosting many ENCODE, ROADMAP, BLUEPRINT, and DEEP samples. These are often hosted as normalized signal tracks that can be used with `bigwigCompare`, `multiBigwigSummary`, and `computeMatrix`. As of version 2.4.0, the aforementioned tools can now access signal files hosted on deepBlue. To do so, simply specify the "experiment name" from deepBlue, such as:

.. code:: bash

    $ bigwigCompare -b1 S002R5H1.ERX300721.H3K4me3.bwa.GRCh38.20150528.bedgraph -b2 S002R5H1.ERX337057.Input.bwa.GRCh38.20150528.bedgraph -p 10 -o bwCompare.bw

The file names given to the aforementioned commands are in the "Name" column in deepBlue. Any file ending in ".wig", ".wiggle", ".bedgraph" or otherwise not present on the file system (and not beginning with "http" or "ftp") is assumed to be hosted on deepBlue. This means that for ENCODE samples, one can simply use the ENCODE ID (e.g., "ENCFF721EKA").

Internally, deepTools queries deepBlue and creates a temporary bigWig file including signal in all of the regions that deepTools will use. By default, these temporary files are deleted after the command finishes. This can be prevented by specifying `--deepBlueKeepTemp`. The directory to which the temporary files are written can be specified by `--deepBlueTempDir`. If you intend to use the same sample multiple times with the same basic command (e.g., computeMatrix with the same regions or bigwigCompare with different samples), then considerable time can be saved by keeping the temporary bigWig file and simply specifying it in subsequent runs (i.e., deepTools won't magically find the previous file, you need to specify it).

Note that some datasets may be restricted access. In such cases, you can request an account and will receive a "user key". You can then provide that to `bigwigCompare`, `multiBigwigSummary`, or `computeMatrix` using the `--userKey` option. In the off-chance that you have access to other deepBlue servers aside from the main one (http://deepblue.mpi-inf.mpg.de/xmlrpc), you can specify that with the `--deepBlueURL` option.

.. warning:: bigwigCompare can be incredibly slow due to essentially downloading entire samples. It's faster to simply download bigWig files from the original source.
