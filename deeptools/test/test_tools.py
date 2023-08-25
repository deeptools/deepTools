from subprocess import PIPE, run
import os
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

TOMLFILE = os.path.dirname(os.path.abspath(__file__)) + "/../../pyproject.toml"


def test_tools():
    """
    Check every script in 'pyproject.toml'
    makes sure the version of all tools == version set in toml file
    makes sure exitcodes are all 0
    """
    with open(TOMLFILE, 'rb') as f:
        _toml = tomllib.load(f)
    for _p in _toml['project']['scripts'].keys():
        _res = run(
            [_p, f"--version"],
            stdout=PIPE,
            stderr=PIPE
        )
        _version = _res.stdout.decode().splitlines()[0]
        assert f"{_version}" == f"{_p} {_toml['project']['version']}"
        assert f"{_res.returncode}" == f"0"

