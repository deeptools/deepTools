import re
import sys
import tomllib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def cargo_rust_version() -> str:
    data = tomllib.loads((REPO_ROOT / "Cargo.toml").read_text())
    return data["package"]["rust-version"]


def pixi_rust_version() -> str:
    data = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())
    spec = data["tool"]["pixi"]["dependencies"]["rust"]
    match = re.fullmatch(r">=\s*([\d.]+)", spec)
    if not match:
        sys.exit(f"pyproject.toml: unexpected [tool.pixi.dependencies] rust spec {spec!r}, "
                  "expected '>=X.Y' so it can be compared to Cargo.toml's rust-version")
    return match.group(1)


def readthedocs_rust_version() -> str:
    text = (REPO_ROOT / ".readthedocs.yaml").read_text()
    match = re.search(r'^\s*rust:\s*"([\d.]+)"', text, re.MULTILINE)
    if not match:
        sys.exit("could not find a 'rust: \"X.Y\"' line under build.tools in .readthedocs.yaml")
    return match.group(1)


def deeptools_version() -> str:
    data = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())
    return data["project"]["version"]


def galaxy_tool_version() -> str:
    text = (REPO_ROOT / "galaxy/wrapper/deepTools_macros.xml").read_text()
    match = re.search(r'<token name="@TOOL_VERSION@">([^<]+)</token>', text)
    if not match:
        sys.exit("could not find a '<token name=\"@TOOL_VERSION@\">' line in "
                  "galaxy/wrapper/deepTools_macros.xml")
    return match.group(1)


def main() -> None:
    cargo_version = cargo_rust_version()
    rust_checks = {
        "pyproject.toml [tool.pixi.dependencies] rust": pixi_rust_version(),
        ".readthedocs.yaml build.tools.rust": readthedocs_rust_version(),
    }
    mismatches = {
        f"rust: {name}": (cargo_version, version)
        for name, version in rust_checks.items()
        if version != cargo_version
    }

    pkg_version = deeptools_version()
    galaxy_version = galaxy_tool_version()
    if galaxy_version != pkg_version:
        mismatches["deeptools: galaxy/wrapper/deepTools_macros.xml @TOOL_VERSION@"] = (
            pkg_version,
            galaxy_version,
        )

    if mismatches:
        print("Version mismatch detected:", file=sys.stderr)
        for name, (expected, actual) in mismatches.items():
            print(f"  {name}: expected {expected!r}, found {actual!r}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
