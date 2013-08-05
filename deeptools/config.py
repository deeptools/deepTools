import os
import ConfigParser
from pkg_resources import resource_stream


"""
If the environment variable DEEP_TOOLS_NO_CONFIG is set the binaries 
from the PATH will be taken. That is used in the Galaxy Tool Shed integration
"""

if  os.environ.get('DEEP_TOOLS_NO_CONFIG', False):
    config = ConfigParser.ConfigParser()
    config.add_section('general')
    config.set('general', 'default_proc_number', 'max/2')

    config.add_section('external_tools')
    config.set('external_tools', 'sort', 'sort')
    config.set('external_tools', 'samtools', 'samtools')
    config.set('external_tools', 'bedgraph_to_bigwig', 'bedGraphToBigWig')
    config.set('external_tools', 'bigwig_info', 'bigWigInfo')

else:
    config = ConfigParser.ConfigParser()
    config_file = resource_stream('config', 'deepTools.cfg')
    config.readfp(config_file)
