import ConfigParser
import os

config = ConfigParser.ConfigParser()
config.read(os.path.dirname(__file__) + "/../config/deepTools.cfg")
