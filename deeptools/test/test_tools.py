import subprocess
import os

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../../bin"


def test_tools():
    """
    Checks everything that is in /bin/
    and tries to run it

    """
    if os.path.isdir(ROOT):
        for _file in os.listdir(ROOT):
            print(_file)
            if os.path.isfile(os.path.join(ROOT, _file)):
                subprocess.check_call("{}/{} --version".format(ROOT, _file).split())
