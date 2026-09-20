import argparse
import gzip
import re
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

DOCS_PATH = Path(__file__).resolve().parent.parent / "docs/content/feature/effectiveGenomeSize.rst"
CACHE_DIR = Path(__file__).resolve().parent / "genomes"
READ_LENGTHS = [50, 75, 100, 150, 200, 250]
REQUIRED_EXECUTABLES = ["faCount", "unique-kmers.py"]

GENOME_URLS = {
    "GRCh37": "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/14/ensembl/2013_09/genome/softmasked.fa.bgz",
    "GRCh38": "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/405/29/ensembl/2026_04/genome/softmasked.fa.bgz",
    "T2T/CHM13CAT_v2": "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/009/914/755/4/ensembl/2022_07/genome/softmasked.fa.bgz",
    "GRCm37 (mm9)": "https://ftp.ensembl.org/pub/release-65/fasta/mus_musculus/dna/Mus_musculus.NCBIM37.65.dna_rm.toplevel.fa.gz",
    "GRCm38 (mm10)": "https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz",
    "GRCm39": "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/635/9/ensembl/2026_04/genome/softmasked.fa.bgz",
    "dm3": "https://ftp.ensembl.org/pub/release-77/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP5.dna_sm.toplevel.fa.gz",
    "dm6": "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/215/4/flybase/2022_07/genome/softmasked.fa.bgz",
    "GRCz10": "https://ftp.ensembl.org/pub/release-91/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna_sm.toplevel.fa.gz",
    "GRCz11": "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/002/035/4/ensembl/2018_04/genome/softmasked.fa.bgz",
    "WBcel235": "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/002/985/3/wormbase/2014_10/genome/softmasked.fa.bgz",
    "TAIR10": "https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/000/001/735/1/community_araport11/2010_09/genome/softmasked.fa.bgz",
}

TABLE_RE = re.compile(r"^\+(?:[-=]+\+)+\n(?:\|.*\|\n\+(?:[-=]+\+)+\n)+", re.MULTILINE)


def url_exists(url: str, timeout: float = 15) -> bool:
    for method, headers in (("HEAD", {}), ("GET", {"Range": "bytes=0-0"})):
        try:
            req = urllib.request.Request(url, method=method, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                if 200 <= resp.status < 400:
                    return True
        except Exception:  # noqa: S110
            pass
    return False


def fetch(name: str, url: str) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    slug = re.sub(r"[^A-Za-z0-9.]+", "_", name)
    dest = CACHE_DIR / f"{slug}__{Path(urllib.parse.urlparse(url).path).name}"
    if not dest.exists():
        print(f"downloading {url} -> {dest}")
        urllib.request.urlretrieve(url, dest)
    if dest.suffix not in (".gz", ".bgz"):
        return dest
    fasta = dest.with_suffix("")
    if not fasta.exists():
        with gzip.open(dest, "rb") as src, open(fasta, "wb") as out:
            shutil.copyfileobj(src, out)
    return fasta


def fa_count(fasta: Path) -> str:
    if shutil.which("faCount") is None:
        raise RuntimeError("faCount not on PATH: http://hgdownload.soe.ucsc.edu/admin/exe/")
    out = subprocess.run(["faCount", str(fasta)], capture_output=True, text=True, check=True).stdout
    total = next(line for line in out.splitlines() if line.startswith("total"))
    _, length, _a, _c, _g, _t, n, *_ = total.split()
    return str(int(length) - int(n))


def unique_kmers(fasta: Path, read_length: int) -> str:
    """`unique-kmers.py -k <read_length> genome.fa` -> estimated unique k-mers."""
    if shutil.which("unique-kmers.py") is None:
        raise RuntimeError("unique-kmers.py not on PATH: pixi run -e ess ...")
    out = subprocess.run(
        ["unique-kmers.py", "-k", str(read_length), str(fasta)],
        capture_output=True, text=True, check=True,
    ).stderr
    match = re.search(r"unique \d+-mers.*?:\s*(\d+)", out, re.IGNORECASE)
    if not match:
        raise RuntimeError(f"couldn't parse unique-kmers.py output:\n{out}")
    return match.group(1)


def parse_grid_table(raw: str):
    rows = [[cell.strip() for cell in line.strip("|").split("|")]
            for line in raw.splitlines() if line.startswith("|")]
    return rows[0], rows[1:]


def rst_grid_table(headers, rows) -> str:
    widths = [max(len(str(cell)) for cell in col) for col in zip(headers, *rows)]

    def rule(char):
        return "+" + "+".join(char * (w + 2) for w in widths) + "+"

    def line(cells):
        return "|" + "|".join(f" {c!s:<{w}} " for c, w in zip(cells, widths)) + "|"

    out = [rule("-"), line(headers), rule("=")]
    for row in rows:
        out += [line(row), rule("-")]
    return "\n".join(out) + "\n"


def process_genome(name: str, url: str):
    fasta = fetch(name, url)
    size = fa_count(fasta)
    kmers = {read_length: unique_kmers(fasta, read_length) for read_length in READ_LENGTHS}
    return name, size, kmers


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--threads", type=int, default=1,
                         help="number of genomes to process concurrently (default: %(default)s)")
    parser.add_argument("--dry-run", action="store_true",
                         help="only check that required executables are on PATH, then exit")
    args = parser.parse_args()

    if args.dry_run:
        missing_exes = [exe for exe in REQUIRED_EXECUTABLES if shutil.which(exe) is None]
        for exe in REQUIRED_EXECUTABLES:
            print(f"{'missing' if exe in missing_exes else 'found  '}: {exe}")

        pending = [(name, url) for name, url in GENOME_URLS.items() if url]
        broken = []
        with ThreadPoolExecutor(max_workers=args.threads) as pool:
            for name, url, ok in pool.map(lambda nu: (*nu, url_exists(nu[1])), pending):
                print(f"{'found  ' if ok else 'missing'}: {name} -> {url}")
                if not ok:
                    broken.append(name)

        sys.exit(1 if missing_exes or broken else 0)

    content = DOCS_PATH.read_text()
    tables = list(TABLE_RE.finditer(content))
    if len(tables) != 2:
        raise RuntimeError(f"expected 2 tables in {DOCS_PATH}, found {len(tables)}")

    _header1, rows1 = parse_grid_table(tables[0].group())
    header2, rows2 = parse_grid_table(tables[1].group())
    sizes = dict(rows1)
    kmer_sizes = {row[0]: dict(zip(header2[1:], row[1:])) for row in rows2}

    pending = [(name, url) for name, url in GENOME_URLS.items() if url]
    try:
        if pending:
            with ThreadPoolExecutor(max_workers=args.threads) as pool:
                for name, size, kmers in pool.map(process_genome, *zip(*pending)):
                    sizes[name] = size
                    for read_length, value in kmers.items():
                        kmer_sizes.setdefault(str(read_length), {})[name] = value
    finally:
        shutil.rmtree(CACHE_DIR, ignore_errors=True)

    names = list(GENOME_URLS)
    table1 = rst_grid_table(["Genome", "Effective size"],
                             [[name, sizes.get(name, "?")] for name in names])
    table2 = rst_grid_table(
        ["Read length", *names],
        [[str(rl), *(kmer_sizes.get(str(rl), {}).get(name, "?") for name in names)]
         for rl in READ_LENGTHS],
    )

    content = (content[:tables[0].start()] + table1 + content[tables[0].end():tables[1].start()]
               + table2 + content[tables[1].end():])
    DOCS_PATH.write_text(content)
    print(f"updated {DOCS_PATH}")


if __name__ == "__main__":
    main()
