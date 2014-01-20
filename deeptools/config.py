import os
import ConfigParser


"""
If the environment variable DEEP_TOOLS_NO_CONFIG is set the binaries
from the PATH will be taken. That is used in the Galaxy Tool Shed integration
"""

if os.environ.get('DEEP_TOOLS_NO_CONFIG', False):
    config = ConfigParser.ConfigParser()
    config.add_section('general')
    config.set('general', 'default_proc_number', 'max/2')

    config.add_section('external_tools')
    config.set('external_tools', 'sort', 'sort')
    config.set('external_tools', 'samtools', 'samtools')
    config.set('external_tools', 'bedgraph_to_bigwig', 'bedGraphToBigWig')
    config.set('external_tools', 'bigwig_info', 'bigWigInfo')

else:
    import pkg_resources

    # load the deepTools configuration file
    # that should be located under the root folder
    # of the deepTools instalation
    config_file = pkg_resources.resource_filename(__name__,
                                                  'config/deeptools.cfg')
    config = ConfigParser.ConfigParser()
    config.readfp(open(config_file, 'r'))
