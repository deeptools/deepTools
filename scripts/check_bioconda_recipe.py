import argparse
import os
import re
import tomllib
import urllib.request
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
RECIPE_URL = "https://raw.githubusercontent.com/bioconda/bioconda-recipes/master/recipes/deeptools/meta.yaml"
PYPI_TO_CONDA = {"matplotlib": "matplotlib-base"}


def load_recipe(recipe_dir: str | None) -> tuple[str, dict]:
    if recipe_dir:
        raw = (Path(recipe_dir) / "meta.yaml").read_text()
    else:
        with urllib.request.urlopen(RECIPE_URL, timeout=15) as response:
            raw = response.read().decode()
    text = re.sub(r"\{%.*?%\}", "", raw)
    text = re.sub(r"\{\{.*?\}\}", "JINJA", text)
    return raw, yaml.safe_load(text)


def split_requirement(requirement: str) -> tuple[str, str]:
    name, spec = re.match(r"([\w.\-]+)\s*(.*)", requirement.strip()).groups()
    spec = re.sub(r"(\.0)+(?=$|,)", "", spec.replace(" ", ""))
    return name.lower().replace("_", "-"), spec


def find_drift(raw: str, recipe: dict, project: dict):
    scripts = project["project"].get("scripts", {})
    build = recipe.get("build") or {}

    entry_points = dict(e.replace(" ", "").split("=", 1) for e in build.get("entry_points") or [])
    for name, target in entry_points.items():
        if name not in scripts:
            yield f"recipe entry point `{name} = {target}` is not in [project.scripts]"
        elif scripts[name] != target:
            yield f"entry point `{name}`: recipe has `{target}`, pyproject.toml `{scripts[name]}`"
    if entry_points:
        for name in sorted(scripts.keys() - entry_points.keys()):
            yield f"`{name}` from [project.scripts] is missing from the recipe's entry_points"

    for command in (recipe.get("test") or {}).get("commands") or []:
        if command.split()[0] not in scripts:
            yield f"recipe test `{command}` is not a deepTools command"

    wanted = {}
    for requirement in project["project"].get("dependencies", []):
        name, spec = split_requirement(requirement)
        wanted[PYPI_TO_CONDA.get(name, name)] = spec
    have = dict(split_requirement(r) for r in recipe["requirements"]["run"])
    have.pop("python", None)
    for name in sorted(wanted.keys() | have.keys()):
        if name not in have:
            yield f"run dependency `{name} {wanted[name]}` is in pyproject.toml but not in the recipe"
        elif name not in wanted:
            yield f"run dependency `{name} {have[name]}` is in the recipe but not in pyproject.toml"
        elif wanted[name] != have[name]:
            yield f"run dependency `{name}`: recipe has `{have[name]}`, pyproject.toml `{wanted[name]}`"

    minimum = re.search(r">=\s*(\d+)\.(\d+)", project["project"]["requires-python"])
    skip = re.search(r"skip:\s*True\s*#\s*\[\s*py\s*<\s*(\d+)\s*\]", raw)
    if minimum and (not skip or skip.group(1) != "".join(minimum.groups())):
        yield (f"requires-python is `{project['project']['requires-python']}`, expected "
               f"`skip: True  # [py < {''.join(minimum.groups())}]` in the recipe")


def main() -> None:
    parser = argparse.ArgumentParser(description="Warn about drift between the bioconda recipe and pyproject.toml.")
    parser.add_argument("--recipe", help="local recipes/deeptools directory (default: bioconda-recipes master)")
    parser.add_argument("--pyproject", default=REPO_ROOT / "pyproject.toml")
    args = parser.parse_args()

    raw, recipe = load_recipe(args.recipe)
    print(f"recipe: {Path(args.recipe) / 'meta.yaml' if args.recipe else RECIPE_URL}")
    project = tomllib.loads(Path(args.pyproject).read_text())
    drift = list(find_drift(raw, recipe, project))

    in_actions = os.environ.get("GITHUB_ACTIONS") == "true"
    for item in drift:
        print(f"::warning title=bioconda recipe drift::{item}" if in_actions else f"drift: {item}")
    if not drift:
        print("no drift between the bioconda recipe and pyproject.toml")

    if summary := os.environ.get("GITHUB_STEP_SUMMARY"):
        lines = [f"- {item}" for item in drift] or ["No drift."]
        with open(summary, "a") as handle:
            handle.write("## bioconda recipe drift\n\n" + "\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
