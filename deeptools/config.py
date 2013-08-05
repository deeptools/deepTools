from os import environ
import ConfigParser
from pkg_resources import resource_stream
config = ConfigParser.ConfigParser()
config_file = resource_stream('config', 'deepTools.cfg')
config.readfp(config_file)

def config_get(class_name, name):
    try:
        # check if an enviromental variable is set.
        # if this is the case, it will override the
        # config file settings
        env_class = class_name + '_env'
        env_name = config.get(env_class, name)
        res = environ[env_name]
    except KeyError:
        res = config.get(class_name, name)

    return res
