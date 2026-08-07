from pathlib import Path
import sys
import tomllib
import matplotlib
import deeptools


def find_pyproject():

    current_path = Path(__file__).resolve()

    for parent in current_path.parents:
        pyproject = parent / "pyproject.toml"

        if pyproject.exists():
            return pyproject

    raise FileNotFoundError(
        "pyproject.toml not found"
    )


def get_matplotlib_requirement():

    pyproject = find_pyproject()

    with open(pyproject, "rb") as f:
        config = tomllib.load(f)

    dependencies = config["project"]["dependencies"]

    for dependency in dependencies:
        if dependency.startswith("matplotlib"):
            return dependency

    raise RuntimeError(
        "matplotlib dependency not found in pyproject.toml"
    )


def test_matplotlib_version():

    print("\nPython executable:")
    print(sys.executable)

    print("\nMatplotlib version:")
    print(matplotlib.__version__)

    print("\nMatplotlib path:")
    print(matplotlib.__file__)

    print("\ndeepTools path:")
    print(deeptools.__file__)

    requirement = get_matplotlib_requirement()

    print("\nMatplotlib requirement from pyproject.toml:")
    print(requirement)

    # Example: "matplotlib >= 3.9"
    min_version = (
        requirement
        .split(">=")[1]
        .strip()
    )

    required = tuple(
        map(int, min_version.split(".")[:2])
    )

    installed = tuple(
        map(int, matplotlib.__version__.split(".")[:2])
    )

    assert installed >= required, (
        f"Installed matplotlib {matplotlib.__version__} "
        f"does not satisfy {requirement}"
    )