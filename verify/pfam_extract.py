#!/usr/bin/env python3
"""Extract the largest families from Pfam-A.seed.gz (Stockholm) as FASTA files and
generate a scale differential manifest.

Stream-parses the gzip, keeps families with at least --min-seqs sequences, writes
the --top largest (by residues = seqs x columns) to verify/corpus/pfam/, and emits
a manifest exercising trimming, statistics and format conversion on each.
"""
import argparse
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def parse_families(path, min_seqs):
    current = None
    name = None
    count = 0
    with gzip.open(path, "rt", encoding="latin1", errors="replace") as handle:
        for line in handle:
            tag = line[:5]
            if line.startswith("//"):
                if current and count >= min_seqs:
                    yield name, current
                current = None
                continue
            if tag == "#=GF ":
                parts = line.split(None, 3)
                if len(parts) >= 3 and parts[1] == "ID":
                    name = parts[2].strip()
            elif tag == "#=GS " or (tag != "#=G" and tag != "#=C" and line[0] != "#"):
                parts = line.split(None, 2)
                if len(parts) >= 2 and not line.startswith("#"):
                    seqname = parts[0]
                    residues = parts[-1].strip().replace(".", "").replace("-", "")
                    if current is None:
                        current = {}
                    if seqname in current:
                        current[seqname] += residues  # wrapped line continuation
                    else:
                        current[seqname] = residues
                    count = len(current)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed-gz", default=str(ROOT / "verify/downloads/Pfam-A.seed.gz"))
    parser.add_argument("--min-seqs", type=int, default=500)
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--out-dir", default=str(ROOT / "verify/corpus/pfam"))
    parser.add_argument("--manifest", default=str(ROOT / "verify/pfam_manifest.txt"))
    options = parser.parse_args()

    out_dir = Path(options.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    families = []
    for name, seqs in parse_families(options.seed_gz, options.min_seqs):
        width = max(len(s) for s in seqs.values())
        cells = len(seqs) * width
        families.append((name, seqs, width, cells))
    families.sort(key=lambda f: -f[3])
    families = families[: options.top]

    rows = []
    for rank, (name, seqs, width, cells) in enumerate(families, 1):
        safe = "".join(c if c.isalnum() else "_" for c in name)[:40]
        path = out_dir / ("pfam_%02d_%s.fasta" % (rank, safe))
        with open(path, "w", encoding="ascii", newline="\n") as handle:
            for seqname, residues in seqs.items():
                handle.write(">%s\n" % seqname)
                for start in range(0, len(residues), 60):
                    handle.write(residues[start:start + 60] + "\n")
        relative = path.relative_to(ROOT).as_posix()
        rows.append((relative, len(seqs), width, cells))

    manifest = ["# Pfam-A.seed scale suite: %d largest families (min %d seqs)"
                % (len(rows), options.min_seqs)]
    for relative, nseq, width, cells in rows:
        manifest.append("-in %s -gappyout -out $OUT" % relative)
        manifest.append("-in %s -automated1 -out $OUT" % relative)
        manifest.append("-in %s -sgc -out $OUT" % relative)
        manifest.append("-in %s -strict -fasta -out $OUT" % relative)
    Path(options.manifest).write_text("\n".join(manifest) + "\n", encoding="utf-8", newline="\n")

    for relative, nseq, width, cells in rows:
        print("%-45s seqs=%-6d cols=%-6d cells=%.1fM" % (relative, nseq, width, cells / 1e6))
    print("manifest cases:", sum(1 for _ in open(options.manifest) if _.strip() and not _.startswith("#")))


if __name__ == "__main__":
    main()
