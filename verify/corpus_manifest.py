#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Record and re-check the bytes of a downloaded-data corpus with SHA-256.

A corpus extracted from Pfam or eggNOG is not committed: it is hundreds of megabytes.
What is committed is its hash list, so that a later extraction can be compared file by
file with the bytes an earlier verification run actually used. The lists are readable
by `sha256sum -c` from the repository root.

    python3 verify/corpus_manifest.py write verify/corpus/eggnog verify/corpus_manifests/eggnog.sha256
    python3 verify/corpus_manifest.py verify verify/corpus_manifests/eggnog.sha256
"""
import argparse
import hashlib
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    sha = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            sha.update(block)
    return sha.hexdigest()


def display(path):
    try:
        return Path(path).resolve().relative_to(ROOT).as_posix()
    except ValueError:  # outside the repository
        return Path(path).resolve().as_posix()


def entries_of(manifest):
    """{path: sha256} for the entries of an existing hash list ({} when absent)."""
    entries = {}
    manifest = Path(manifest)
    if manifest.exists():
        for line in manifest.read_text(encoding="utf-8").splitlines():
            if line.strip() and not line.startswith("#"):
                expected, relative = line.split(None, 1)
                entries[relative.strip()] = expected
    return entries


def differences_between(manifest, new):
    """Differences between an existing hash list and a {path: sha256} mapping of the
    bytes an extraction is about to write. Paths of the list that the extraction does
    not produce are ignored (other files kept in the corpus)."""
    old = entries_of(manifest)
    differences = []
    for relative in sorted(new):
        if relative not in old:
            differences.append("NEW " + relative)
        elif old[relative] != new[relative]:
            differences.append("DIFFERS " + relative)
    return differences


def write_hashes(directory, manifest, notes=()):
    directory, manifest = Path(directory), Path(manifest)
    files = sorted(p for p in directory.rglob("*") if p.is_file())
    if not files:
        raise SystemExit("no files under " + str(directory))
    manifest.parent.mkdir(parents=True, exist_ok=True)
    lines = ["# corpus: %s" % display(directory)]
    lines += ["# %s" % note for note in notes]
    lines += ["# %d file(s), SHA-256, path relative to the repository root" % len(files)]
    lines += ["%s  %s" % (digest(p), display(p)) for p in files]
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    total = sum(p.stat().st_size for p in files)
    print("HASHES %s files=%d bytes=%d" % (display(manifest), len(files), total))
    return manifest


def verify_hashes(manifest):
    manifest = Path(manifest)
    missing, differing, checked = [], [], 0
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        expected, relative = line.split(None, 1)
        path = ROOT / relative.strip()
        checked += 1
        if not path.is_file():
            missing.append(relative.strip())
        elif digest(path) != expected:
            differing.append(relative.strip())
    for relative in missing:
        print("MISSING %s" % relative)
    for relative in differing:
        print("DIFFERS %s" % relative)
    print("VERIFY %s checked=%d missing=%d differing=%d" % (
        display(manifest), checked, len(missing), len(differing)))
    return 1 if missing or differing else 0


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    write = sub.add_parser("write")
    write.add_argument("directory", type=Path)
    write.add_argument("manifest", type=Path)
    write.add_argument("--note", action="append", default=[])
    check = sub.add_parser("verify")
    check.add_argument("manifest", type=Path)
    options = parser.parse_args()
    if options.command == "write":
        write_hashes(options.directory, options.manifest, options.note)
        return 0
    return verify_hashes(options.manifest)


if __name__ == "__main__":
    sys.exit(main())
