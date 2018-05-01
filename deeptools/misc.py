import os

# This should force numpy to run single threaded. See issue #697
# This module MUST be imported before numpy
# Note that these environment variables are internal to deepTools (they won't exist on the shell after the command completes)
if not 'MKL_NUM_THREADS' in os.environ:
    os.environ['MKL_NUM_THREADS'] = 'sequential'
if not 'NUMEXPR_NUM_THREADS' in os.environ:
    os.environ['NUMEXPR_NUM_THREADS'] = '1'
if not 'OMP_NUM_THREADS' in os.environ:
    os.environ['OMP_NUM_THREADS'] = '1'
if not 'VECLIB_MAXIMUM_THREADS' in os.environ:
    os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
