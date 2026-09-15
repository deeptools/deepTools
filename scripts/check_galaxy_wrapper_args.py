import ast
import re
import sys
import tomllib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
WRAPPER_DIR = REPO_ROOT / "galaxy" / "wrapper"
DEEPTOOLS_SRC = REPO_ROOT / "pydeeptools" / "deeptools"
MACROS_FILE = WRAPPER_DIR / "deepTools_macros.xml"
PYPROJECT_FILE = REPO_ROOT / "pyproject.toml"

TOKEN_RE = re.compile(r'<token name="(@[A-Za-z0-9_]+@)">(.*?)</token>', re.S)
BINARY_TOKEN_RE = re.compile(r'<token name="@BINARY@">([^<]+)</token>')
COMMAND_RE = re.compile(
    r"<command[^>]*>\s*(?:<!\[CDATA\[)?(.*?)(?:\]\]>)?\s*</command>", re.S
)
REAL_ATTR_FLAG_RE = re.compile(r'(?:truevalue|falsevalue)="(--?[A-Za-z][\w-]*)"')
DOC_ATTR_FLAG_RE = re.compile(r'argument="(--?[A-Za-z][\w-]*)"')
FLAG_TOKEN_RE = re.compile(
    r"(?<![\w-])--[A-Za-z][\w-]*|(?<![\w-])-[A-Za-z][A-Za-z0-9]*"
)


def load_tokens(text):
    return {name: value for name, value in TOKEN_RE.findall(text)}


def expand_tokens(text, tokens):
    for _ in range(5):
        new_text = text
        for name, value in tokens.items():
            new_text = new_text.replace(name, value)
        if new_text == text:
            break
        text = new_text
    return text


def extract_command_flags(command_text, binary):
    """
    Ignore flags belonging to `ln`, `samtools`, etc.
    """
    flags = set()
    for segment in command_text.split("&&"):
        seg = segment.strip()
        words = seg.split(None, 1)
        if not words or words[0] != binary:
            continue
        flags |= set(FLAG_TOKEN_RE.findall(seg))
    return flags


def wrapper_flags(xml_path, macro_tokens):
    """Returns (binary, real_flags, doc_only_flags)."""
    text = xml_path.read_text()
    binary_match = BINARY_TOKEN_RE.search(text)
    if not binary_match:
        return None, set(), set()
    binary = binary_match.group(1).strip()

    tokens = dict(macro_tokens)
    tokens.update(load_tokens(text))

    command_match = COMMAND_RE.search(text)
    command_text = expand_tokens(command_match.group(1), tokens) if command_match else ""

    real_flags = extract_command_flags(command_text, binary)
    real_flags |= {m.group(1) for m in REAL_ATTR_FLAG_RE.finditer(text)}

    doc_flags = {m.group(1) for m in DOC_ATTR_FLAG_RE.finditer(text)} - real_flags
    return binary, real_flags, doc_flags


def load_entry_points():
    data = tomllib.loads(PYPROJECT_FILE.read_text())
    scripts = data["project"]["scripts"]
    return {name: target.split(":", 1)[0] for name, target in scripts.items()}


def module_path_for(dotted):
    parts = dotted.split(".")
    if parts[0] != "deeptools" or len(parts) < 2:
        return None
    return DEEPTOOLS_SRC / Path(*parts[1:]).with_suffix(".py")


def local_deeptools_imports(tree):
    """
    python imports
    """
    mods = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("deeptools"):
            for alias in node.names:
                candidate = f"{node.module}.{alias.name}"
                if module_path_for(candidate) and module_path_for(candidate).is_file():
                    mods.add(candidate)
                else:
                    mods.add(node.module)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith("deeptools") and module_path_for(alias.name):
                    mods.add(alias.name)
    return mods


def collect_flags(dotted_module, visited=None):
    """Statically walk add_argument(...) calls reachable from a module,
    following its local `deeptools.*` imports (e.g. parserCommon)."""
    if visited is None:
        visited = set()
    if dotted_module in visited:
        return set()
    visited.add(dotted_module)

    path = module_path_for(dotted_module)
    if path is None or not path.is_file():
        return set()

    tree = ast.parse(path.read_text(), filename=str(path))
    flags = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "add_argument":
            flags.update(a.value for a in node.args if isinstance(a, ast.Constant) and isinstance(a.value, str))

    for mod in local_deeptools_imports(tree):
        flags |= collect_flags(mod, visited)
    return flags


def main():
    entry_points = load_entry_points()
    macro_tokens = load_tokens(MACROS_FILE.read_text())

    errors = []
    doc_warnings = []

    for xml_path in sorted(WRAPPER_DIR.glob("*.xml")):
        if xml_path == MACROS_FILE:
            continue

        binary, real_flags, doc_flags = wrapper_flags(xml_path, macro_tokens)
        if binary is None:
            continue

        module = entry_points.get(binary)
        if module is None:
            errors.append(f"{xml_path.name}: no [project.scripts] entry point for binary '{binary}'")
            continue

        valid_flags = collect_flags(module)
        if not valid_flags:
            errors.append(f"{xml_path.name}: could not statically find any add_argument() for '{binary}' ({module})")
            continue

        for flag in sorted(real_flags):
            if flag not in valid_flags:
                errors.append(f"{xml_path.name}: '{flag}' is not a recognized argument of {binary} ({module})")

        for flag in sorted(doc_flags):
            if flag not in valid_flags:
                doc_warnings.append(
                    f"{xml_path.name}: argument=\"{flag}\" doesn't match any flag of {binary} "
                    f"({module}) -- doc-only, not used on the command line"
                )

    if doc_warnings:
        print("Stale argument= attributes (documentation only, not passed on the command line):")
        for w in doc_warnings:
            print(f"  - {w}")
        print()

    if errors:
        print("Argument mismatches between galaxy/wrapper/*.xml and the Python CLI:")
        for e in errors:
            print(f"  - {e}")
        print(f"\n{len(errors)} error(s) found.")
        return 1

    print("All galaxy wrapper arguments match their Python CLI definitions.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
