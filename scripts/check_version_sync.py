import json
import re
import sys
import tomllib
import urllib.error
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# RTD only ships specific Rust toolchains, so build.tools.rust can't be set to an
# arbitrary version string the way Cargo.toml's rust-version can. The allowed set
# is published as the JSON schema RTD uses to validate .readthedocs.yaml.
READTHEDOCS_SCHEMA_URL = (
    "https://raw.githubusercontent.com/readthedocs/readthedocs.org/main/"
    "readthedocs/rtd_tests/fixtures/spec/v2/schema.json"
)


def _version_tuple(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


def _find_rust_enum(node: object) -> list[str] | None:
    """Walk the RTD JSON schema looking for the `rust` tool's `enum` list."""
    if isinstance(node, dict):
        rust = node.get("rust")
        if isinstance(rust, dict) and isinstance(rust.get("enum"), list):
            return rust["enum"]
        for value in node.values():
            found = _find_rust_enum(value)
            if found is not None:
                return found
    elif isinstance(node, list):
        for item in node:
            found = _find_rust_enum(item)
            if found is not None:
                return found
    return None


def fetch_readthedocs_allowed_rust_versions() -> list[str] | None:
    """Return RTD's allowed build.tools.rust values, or None if unreachable/unparseable."""
    try:
        with urllib.request.urlopen(READTHEDOCS_SCHEMA_URL, timeout=5) as response:
            schema = json.load(response)
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        print(
            f"warning: could not fetch RTD's schema ({READTHEDOCS_SCHEMA_URL}) to validate "
            f"build.tools.rust ({exc}); skipping the allowed-toolchain check",
            file=sys.stderr,
        )
        return None
    versions = _find_rust_enum(schema)
    if versions is None:
        print(
            f"warning: RTD's schema ({READTHEDOCS_SCHEMA_URL}) no longer has the expected "
            "build.tools.rust enum shape; skipping the allowed-toolchain check",
            file=sys.stderr,
        )
    return versions


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
    match = re.search(r'^\s*rust:\s*"([\w.]+)"', text, re.MULTILINE)
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
    mismatches: dict[str, tuple[str, str]] = {}

    pixi_version = pixi_rust_version()
    if pixi_version != cargo_version:
        mismatches["rust: pyproject.toml [tool.pixi.dependencies] rust"] = (cargo_version, pixi_version)

    rtd_version = readthedocs_rust_version()
    allowed_rtd_versions = fetch_readthedocs_allowed_rust_versions()
    if allowed_rtd_versions is not None and rtd_version not in allowed_rtd_versions:
        mismatches["rust: .readthedocs.yaml build.tools.rust (not an RTD-provided toolchain)"] = (
            f"one of {allowed_rtd_versions}",
            rtd_version,
        )
    elif rtd_version != "latest" and _version_tuple(rtd_version) < _version_tuple(cargo_version):
        mismatches["rust: .readthedocs.yaml build.tools.rust"] = (f">= {cargo_version}", rtd_version)

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
