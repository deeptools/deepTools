Accessing datasets hosted on deepBlue
=====================================

`deepBlue <http://dx.doi.org/10.1093/nar/gkw211>` is an epigenome dataset server hosting many ENCODE, ROADMAP, BLUEPRINT, and DEEP samples. These are often hosted as normalized signal tracks that can be used with `bigwigCompare`, `multiBigwigSummary`, and `computeMatrix`. As of version 2.4.0, the aforementioned tools can now access signal files hosted on deepBlue. To do so, simply specify the "experiment name" from deepBlue, such as:

.. code:: bash

    $ bigwigCompare -b1 S002R5H1.ERX300721.H3K4me3.bwa.GRCh38.20150528.bedgraph -b2 S002R5H1.ERX337057.Input.bwa.GRCh38.20150528.bedgraph -p 10 -o bwCompare.bw

For `multiBigwigSummary`, performance is terrible unless the bin size is increased. Using `-bs 1000000` is suggested in this case. For computeMatrix, it is suggested to remove any unneccessary regions. In general, performance is poor whenever the regions of interest are small (below ~100kb).

Note that some datasets may be restricted access. In such cases, you can request an account and will receive a "user key". You can then provide that to `bigwigCompare`, `multiBigwigSummary`, or `computeMatrix` using the `--userKey` option. In the off-chance that you have access to other deepBlue servers aside from the main one (http://deepblue.mpi-inf.mpg.de/xmlrpc), you can specify that with the `--deepBlueURL` option.
