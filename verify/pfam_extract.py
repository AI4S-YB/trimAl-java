#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Extract the largest families from Pfam-A.seed.gz (Stockholm) as aligned FASTA files
and generate a scale differential manifest.

Stream-parses the gzip, keeps families with at least --min-seqs sequences, takes
the --top largest (by residues = seqs x columns) and emits a manifest exercising
trimming, statistics and format conversion on each.

The alignment is preserved. Pfam marks columns that are insertions relative to the
model with `.` in every sequence that lacks the inserted residues; such columns are
dropped, and the `-` gaps of the remaining columns are kept, so trimAl sees an
aligned family. (An earlier version stripped every gap character, which produced
unaligned sequences: trimAl rejected all of them with "The sequences in the input
alignment should be aligned", and the whole gate compared that one error message.)

The FASTA bytes are compared with the committed hash list
(verify/corpus_manifests/pfam.sha256) before anything is written. A corpus that
differs from the one earlier verification runs used is reported and not written
unless --update-manifest is given, so a changed upstream release cannot silently
replace the evidence.
"""
import argparse
import gzip
import hashlib
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "verify"))

import corpus_manifest  # noqa: E402

HASH_MANIFEST = ROOT / "verify/corpus_manifests/pfam.sha256"


def parse_families(path, min_seqs):
    """(name, {sequence: aligned residues with '.' and '-' kept}) per family."""
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
                    residues = parts[-1].strip()
                    if current is None:
                        current = {}
                    if seqname in current:
                        current[seqname] += residues  # wrapped line continuation
                    else:
                        current[seqname] = residues
                    count = len(current)


def aligned(seqs):
    """Drop the insertion columns (a '.' in any sequence); keep the '-' gaps."""
    width = min(len(seq) for seq in seqs.values())
    if len({len(seq) for seq in seqs.values()}) != 1:
        return None
    keep = [column for column in range(width)
            if all(seq[column] != "." for seq in seqs.values())]
    return {sequence: "".join(seq[column] for column in keep) for sequence, seq in seqs.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed-gz", default=str(ROOT / "verify/downloads/Pfam-A.seed.gz"))
    parser.add_argument("--min-seqs", type=int, default=500)
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--out-dir", default=str(ROOT / "verify/corpus/pfam"))
    parser.add_argument("--manifest", default=str(ROOT / "verify/pfam_manifest.txt"))
    parser.add_argument("--update-manifest", action="store_true",
                        help="accept family bytes that differ from the committed hash list")
    options = parser.parse_args()

    out_dir = Path(options.out_dir)

    # (cells, order, name, seqs, width); pruned to the leaders so that only a few
    # families are ever held in memory
    kept = []
    for order, (name, seqs) in enumerate(parse_families(options.seed_gz, options.min_seqs)):
        seqs = aligned(seqs)
        if seqs is None:
            print("skipped family with unequal sequence lengths:", name)
            continue
        width = len(next(iter(seqs.values())))
        kept.append((len(seqs) * width, order, name, seqs, width))
        if len(kept) > 4 * options.top:
            kept.sort(key=lambda f: (-f[0], f[1]))
            del kept[options.top:]
    kept.sort(key=lambda f: (-f[0], f[1]))
    families = kept[: options.top]

    rows = []
    texts = {}
    for rank, (cells, _, name, seqs, width) in enumerate(families, 1):
        safe = "".join(c if c.isalnum() else "_" for c in name)[:40]
        path = out_dir / ("pfam_%02d_%s.fasta" % (rank, safe))
        text = "".join(">%s\n%s" % (seqname, "".join(
            residues[start:start + 60] + "\n" for start in range(0, len(residues), 60)))
            for seqname, residues in seqs.items())
        relative = path.relative_to(ROOT).as_posix()
        texts[relative] = text
        rows.append((relative, len(seqs), width, cells))

    differences = corpus_manifest.differences_between(
        HASH_MANIFEST, {relative: hashlib.sha256(text.encode("ascii")).hexdigest()
                        for relative, text in texts.items()})
    for line in differences:
        print("CORPUS " + line)
    if differences and HASH_MANIFEST.exists() and not options.update_manifest:
        print("the extracted families differ from the committed hash list; nothing was "
              "written. Re-run with --update-manifest to accept the new bytes",
              file=sys.stderr)
        return 1

    out_dir.mkdir(parents=True, exist_ok=True)
    for relative, text in texts.items():
        (ROOT / relative).write_bytes(text.encode("ascii"))
    corpus_manifest.write_hashes(out_dir, HASH_MANIFEST, [
        "corpus of the Pfam-A.seed gate",
        "source: %s; %d largest families with at least %d sequences"
        % (corpus_manifest.display(options.seed_gz), len(rows), options.min_seqs)])

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
    return 0


if __name__ == "__main__":
    sys.exit(main())
