#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Extract the mega-families of a Pfam-A.full window as FASTA files and generate the
manifest of the mega-family scale gate.

The window is a byte-prefix of Pfam-A.full.gz, so the gzip stream ends mid-family;
families completed before the cut are kept (EOFError tolerated). Every family with at
least --min-seqs sequences is written to a staging directory as soon as it completes
(only one family is ever held in memory); the --top largest by sequence count are
then moved into verify/corpus/pffull/.

Only O(n x L) methods are emitted for the manifest: gappyout / noallgaps / format
conversion / gap statistics. (strict/automated1 need the O(n^2) identity matrix,
which is infeasible at these sizes in either implementation.)

Before a family replaces a corpus file its SHA-256 is compared with the committed
hash list (verify/corpus_manifests/pffull.sha256). Bytes that differ from the ones
earlier verification runs used are reported and discarded unless --update-manifest
is given.
"""
import argparse
import gzip
import hashlib
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "verify"))

import corpus_manifest  # noqa: E402

HASH_MANIFEST = ROOT / "verify/corpus_manifests/pffull.sha256"
METHODS = ["-gappyout", "-noallgaps", "-fasta", "-sgc", "-sgt"]


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


def write_fasta(path, seqs):
    with open(path, "w", encoding="ascii", newline="\n") as handle:
        for seqname, residues in seqs.items():
            handle.write(">%s\n" % seqname)
            for start in range(0, len(residues), 60):
                handle.write(residues[start:start + 60] + "\n")


def digest(path):
    sha = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            sha.update(block)
    return sha.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--window-gz", default=str(ROOT / "verify/downloads/pffull_window.gz"))
    parser.add_argument("--min-seqs", type=int, default=15000)
    parser.add_argument("--top", type=int, default=6)
    parser.add_argument("--out-dir", default=str(ROOT / "verify/corpus/pffull"))
    parser.add_argument("--manifest", default=str(ROOT / "verify/pffull_manifest.txt"))
    parser.add_argument("--update-manifest", action="store_true",
                        help="accept family bytes that differ from the committed hash list")
    options = parser.parse_args()

    out_dir = Path(options.out_dir)
    staging = out_dir / ".staging"
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)

    def safe(name):
        return "".join(c if c.isalnum() else "_" for c in name)[:40]

    found = []
    for name, seqs in iter_families(options.window_gz):
        if len(seqs) < options.min_seqs:
            continue
        width = max(len(s) for s in seqs.values())
        staged = staging / ("%05d_%s.fasta" % (len(found), safe(name)))
        write_fasta(staged, seqs)
        found.append((name, staged, len(seqs), width))
        print("family %-22s seqs=%-7d cols=%-5d cells=%.1fM" % (
            name, len(seqs), width, len(seqs) * width / 1e6), flush=True)
        del seqs
    found.sort(key=lambda f: -f[2])
    top = found[: options.top]
    if not top:
        shutil.rmtree(staging)
        raise SystemExit("no family with >= %d sequences found in window" % options.min_seqs)

    rows = []
    for rank, (name, staged, nseq, width) in enumerate(top, 1):
        final = out_dir / ("pffull_%02d_%s.fasta" % (rank, safe(name)))
        rows.append((final, staged, nseq, width, digest(staged)))

    differences = corpus_manifest.differences_between(
        HASH_MANIFEST, {final.relative_to(ROOT).as_posix(): sha for final, _, _, _, sha in rows})
    for line in differences:
        print("CORPUS " + line)
    if differences and HASH_MANIFEST.exists() and not options.update_manifest:
        shutil.rmtree(staging)
        print("the extracted families differ from the committed hash list; the corpus "
              "was left untouched. Re-run with --update-manifest to accept the new bytes",
              file=sys.stderr)
        return 1

    out_dir.mkdir(parents=True, exist_ok=True)
    for final, staged, _, _, _ in rows:
        staged.replace(final)
    shutil.rmtree(staging)
    corpus_manifest.write_hashes(out_dir, HASH_MANIFEST, [
        "corpus of the Pfam-A.full mega-family gate",
        "source: %s; %d largest families with at least %d sequences"
        % (corpus_manifest.display(options.window_gz), len(rows), options.min_seqs)])

    lines = ["# Pfam-A.full mega families: O(n*L) methods only (identity matrix infeasible at this scale)"]
    for final, _, _, _, _ in rows:
        lines.extend("-in %s %s -out $OUT" % (final.relative_to(ROOT).as_posix(), method)
                     for method in METHODS)
    Path(options.manifest).write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")

    for final, _, nseq, width, _ in rows:
        print("%-50s seqs=%-7d cols=%-5d cells=%.1fM" % (
            final.relative_to(ROOT).as_posix(), nseq, width, nseq * width / 1e6))
    print("manifest cases:", len(lines) - 1)
    print("HASH LIST " + HASH_MANIFEST.relative_to(ROOT).as_posix())
    return 0


if __name__ == "__main__":
    sys.exit(main())
