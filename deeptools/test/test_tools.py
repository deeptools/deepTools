import subprocess
import os

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../../bin"


def test_tools():
    """
    Checks everything that is in /bin/
    and tries to run it

    """
    for file in os.listdir(ROOT):
        print file
        if os.path.isfile(os.path.join(ROOT, file)):
            subprocess.check_call("{}/{} --version".format(ROOT, file).split())


