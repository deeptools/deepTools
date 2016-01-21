bamCompare
===========

.. argparse::
   :ref: deeptools.bamCompare.parseArguments
   :prog: bamCompare

.. warning:: The filtering by deepTools is done *after* the scaling factors are calculated!

    If you know that your files will be strongly affected by the kind of filtering you would like to apply (e.g., removal of duplicates with ``--ignoreDuplicates`` or ignoring reads of low quality) then consider removing those reads *beforehand*. 
