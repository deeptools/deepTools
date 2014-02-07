import os
import sys
import ConfigParser
import subprocess


def checkProgram(program, args, where_to_download):
    """
    deeptools relies on some command line programs
    to work properly. This is a generic routine
    that will check for such programs

    """
    if os.environ.get('DEEP_TOOLS_NO_CONFIG', False):
        return

    try:
        _out = subprocess.Popen([program, args], stderr=subprocess.PIPE,
                                stdout=subprocess.PIPE)
        return True
    except EnvironmentError:
        # handle file not found error.
        # the config file is installed in:
        msg = "\n######################################################\n"\
              "\nThe program *{}* was not found in your PATH. In\n" \
              "order for deeptools to work properly this program needs\n"\
              "to be installed. If you already have a copy of this\n"\
              "program please be sure that it is found in your PATH or\n"\
              "that is referred in the configuration file of deepTools\n"\
              "located at:\n\n{}\n\n" \
              "The program can be downloaded from here:\n " \
              " {}\n\n" \
              "\n########################################################"\
              "\n\n".format(program, config_file, where_to_download)
        sys.stderr.write(msg)

    except Exception as e:
        sys.stderr.write("Error: {}".format(e))

    return False


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

