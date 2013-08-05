import ConfigParser
from pkg_resources import resource_stream
config = ConfigParser.ConfigParser()
config_file = resource_stream('config', 'deepTools.cfg')
config.readfp(config_file)
