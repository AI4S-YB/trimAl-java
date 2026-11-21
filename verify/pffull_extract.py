#!/usr/bin/env python3
"""Extract mega-families (tens of thousands of sequences) from a Pfam-A.full window.

The window is a byte-prefix of Pfam-A.full.gz, so the gzip stream ends mid-family;
families completed before the cut are kept (EOFError tolerated). Only O(n x L)
methods are emitted for the manifest: gappyout / noallgaps / format conversion /
gap statistics. (strict/automated1 need the O(n^2) identity matrix, which is
infeasible at these sizes in either implementation.)
"""
import argparse
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def iter_families(path):
    name = None
    seqs = {}
    with gzip.open(path, "rt", encoding="latin1", errors="replace") as handle:
        try:
            for line in handle:
                if line.startswith("//"):
                    if seqs:
                        yield name, seqs
                    name, seqs = None, {}
                    continue
                if line.startswith("#"):
                    if line.startswith("#=GF"):
                        parts = line.split(None, 3)
                        if len(parts) >= 3 and parts[1] == "ID":
                            name = parts[2].strip()
                    continue
                parts = line.split(None, 2)
                if len(parts) >= 2:
                    seqname, residues = parts[0], parts[-1].strip()
                    seqs[seqname] = seqs.get(seqname, "") + residues
        except (EOFError, OSError):
            pass  # truncated tail of the prefix window


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--window-gz", default=str(ROOT / "verify/downloads/pffull_window.gz"))
    parser.add_argument("--min-seqs", type=int, default=15000)
    parser.add_argument("--top", type=int, default=6)
    parser.add_argument("--out-dir", default=str(ROOT / "verify/corpus/pffull"))
    parser.add_argument("--manifest", default=str(ROOT / "verify/pffull_manifest.txt"))
    options = parser.parse_args()

    out_dir = Path(options.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    families = []
    for name, seqs in iter_families(options.window_gz):
        if len(seqs) >= options.min_seqs:
            width = max(len(s) for s in seqs.values())
            families.append((name, seqs, width))
    families.sort(key=lambda f: -len(f[1]))
    families = families[: options.top]
    if not families:
        raise SystemExit("no family with >= %d sequences found in window" % options.min_seqs)

    rows = []
    for rank, (name, seqs, width) in enumerate(families, 1):
        safe = "".join(c if c.isalnum() else "_" for c in name)[:40]
        path = out_dir / ("pffull_%02d_%s.fasta" % (rank, safe))
        with open(path, "w", encoding="ascii", newline="\n") as handle:
            for seqname, residues in seqs.items():
                handle.write(">%s\n" % seqname)
                for start in range(0, len(residues), 60):
                    handle.write(residues[start:start + 60] + "\n")
        rows.append((path.relative_to(ROOT).as_posix(), len(seqs), width))

    lines = ["# Pfam-A.full mega families: O(n*L) methods only (identity matrix infeasible at this scale)"]
    for relative, _, _ in rows:
        lines.append("-in %s -gappyout -out $OUT" % relative)
        lines.append("-in %s -noallgaps -out $OUT" % relative)
        lines.append("-in %s -fasta -out $OUT" % relative)
        lines.append("-in %s -sgc -out $OUT" % relative)
        lines.append("-in %s -sgt -out $OUT" % relative)
    Path(options.manifest).write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")

    for relative, nseq, width in rows:
        print("%-50s seqs=%-7d cols=%-5d cells=%.1fM" % (relative, nseq, width, nseq * width / 1e6))
    print("manifest cases:", len(lines) - 1)


if __name__ == "__main__":
    main()
