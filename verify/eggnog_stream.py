#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Extract the pinned member alignments from the eggNOG5 Bacteria trimmed-alignment tar.

Reads the archive as one byte stream, either through sequential HTTP Range requests
against the remote tar or from a local prefix of it (`--from-file`), walks its
512-byte tar headers in order and saves the requested complete members. Data of
members that are not wanted is read through when it is short and skipped by resuming
the stream after it when it is long, so only headers and selected members are kept.
A member cut by the end of the available data is left out: its data would be
incomplete.

By default the members are the pinned list below. The corpus was first extracted by
an earlier streaming prototype whose chunk-boundary handling made the selection depend
on download timing, so the eight members are pinned by name to make the corpus
reproducible. `--by-size` selects by stored size instead (for exploring).

Each saved member is decompressed (.faa.gz -> .faa). The corpus is then compared
with the committed hash list (verify/corpus_manifests/eggnog.sha256): a corpus that
differs from the one earlier verification runs used is reported, and the list is not
rewritten unless --update-manifest is given.

    python3 verify/eggnog_stream.py                                  # remote, HTTP Range
    python3 verify/eggnog_stream.py --from-file verify/downloads/eggnog_window.bin
    python3 verify/eggnog_stream.py --verify                         # corpus vs. hash list
    bash verify/pdl.sh <tar-url> <bytes> <out> 16 <workdir>          # ranged download of a prefix

The same extractor drives the 32-case gate manifest verify/eggnog_manifest.txt.
"""
import argparse
import gzip
import hashlib
from pathlib import Path
import sys
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "verify"))

import corpus_manifest  # noqa: E402

URL = "http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2/2_trimmed_algs.tar"
CORPUS = ROOT / "verify/corpus/eggnog"
HASH_MANIFEST = ROOT / "verify/corpus_manifests/eggnog.sha256"
CASE_MANIFEST = ROOT / "verify/eggnog_manifest.txt"
CHUNK = 16 * 1024 * 1024
MIB = 1024 * 1024
METHODS = ["-gappyout", "-automated1", "-strict -fasta", "-sgc"]
# The corpus of every eggNOG gate run so far (see corpus_manifests/eggnog.sha256).
PINNED = ["COG0006", "COG0008", "COG0013", "COG0018", "COG0019", "COG0025", "COG0031",
          "COG0037"]


class RemoteSource:
    """Byte ranges of the remote tar, fetched with HTTP Range requests."""

    size = None

    def __init__(self, url):
        self.url = url

    def read(self, start, end):
        last = None
        for attempt in range(5):
            try:
                request = urllib.request.Request(
                    self.url, headers={"Range": "bytes=%d-%d" % (start, end)})
                with urllib.request.urlopen(request, timeout=180) as response:
                    return response.read()
            except Exception as error:  # transient proxy / network errors
                last = error
                print("  fetch retry %d after error: %r" % (attempt + 1, error),
                      file=sys.stderr, flush=True)
                time.sleep(2 * (attempt + 1))
        raise RuntimeError("range fetch failed: %r" % last)


class LocalSource:
    """Byte ranges of a local file: a prefix window of the tar or a full download."""

    def __init__(self, path):
        self.handle = open(path, "rb")
        self.handle.seek(0, 2)
        self.size = self.handle.tell()

    def read(self, start, end):
        self.handle.seek(start)
        return self.handle.read(end - start + 1)


def octal(raw):
    text = raw.rstrip(b"\0 ").decode("latin1")
    return int(text, 8) if text else 0


def iter_members(source, stop_at, select):
    """Yield (name, size, data) for the complete tar members before offset stop_at
    for which select(name, size) is true.

    Every header is read; the data of a member that is not selected is skipped by
    resuming the stream after it, which costs nothing on a range-readable source.
    `base` is the archive offset of buf[0] and `pos` the index of the next unread
    byte inside buf, so the buffer is compacted only when it runs short.
    """
    buf, pos, base = b"", 0, 0
    long_name = None

    def fill(needed):
        nonlocal buf, pos, base
        if pos and len(buf) - pos < needed:
            buf, base, pos = buf[pos:], base + pos, 0
        while len(buf) < needed:
            start = base + len(buf)
            if start >= stop_at:
                return False
            chunk = source.read(start, min(start + CHUNK, stop_at) - 1)
            if not chunk:
                return False
            buf += chunk
        return True

    while True:
        if not fill(512):
            return
        header = buf[pos:pos + 512]
        pos += 512
        if not header.strip(b"\0"):
            continue              # end-of-archive zero block or padding
        size = octal(header[124:136])
        padded = ((size + 511) // 512) * 512
        kind = header[156:157]
        if kind in (b"L", b"x", b"g", b"K"):
            if not fill(padded):
                return
            data, pos = buf[pos:pos + size], pos + padded
            if kind == b"L":      # GNU long name: belongs to the next header
                long_name = data.rstrip(b"\0")
            continue              # pax / global metadata, not an alignment
        name = (long_name or header[0:100].split(b"\0")[0]).decode("latin1")
        long_name = None
        if not select(name, size):
            end = pos + padded
            if end <= len(buf):
                pos = end                  # already buffered: just advance
            elif padded <= 2 * CHUNK:
                if not fill(padded):       # small member: reading through is cheaper
                    return                 # than a new range request
                pos += padded              # fill() may compact the buffer: advance relatively
            else:                          # large member: resume the stream after it
                buf, base, pos = b"", base + end, 0
            continue
        if not fill(padded):
            return                # the available data ends inside this member
        data = buf[pos:pos + size]
        pos += padded
        yield name, size, data


def extract(source, options):
    stop_at = options.stop_after_mb * MIB
    if source.size is not None:
        stop_at = min(stop_at, source.size)
    lo, hi = options.min_mb * MIB, options.max_mb * MIB
    pinned = None if options.by_size else set(options.members.split(","))
    found, seen, skipped = set(), 0, 0

    def select(name, size):
        nonlocal seen, skipped
        if pinned is not None:
            return Path(name).name.split(".")[0] in pinned - found
        if not lo <= size <= hi:
            return False
        seen += 1
        if len(found) >= options.want:
            return False
        if seen <= options.skip:
            skipped += 1
            return False
        return True

    saved = {}
    for name, size, data in iter_members(source, stop_at, select):
        plain = gzip.decompress(data)
        base = "verify/corpus/eggnog/" + Path(name).name
        saved[base] = (data, hashlib.sha256(data).hexdigest())
        saved[base[:-3]] = (plain, hashlib.sha256(plain).hexdigest())
        found.add(Path(name).name.split(".")[0])
        print("extracted %-38s %7.2f MB gz -> %8.2f MB faa" % (
            Path(name).name, size / MIB, len(plain) / MIB), flush=True)
    print("SUMMARY %s extracted=%d files=%d" % (
        "members_in_window=%d skipped=%d" % (seen, skipped) if pinned is None
        else "pinned=%d" % len(pinned), len(found), len(saved)), flush=True)
    if pinned is not None and pinned - found:
        raise SystemExit("members not found before the stop offset: " + ", ".join(sorted(pinned - found)))
    if not saved:
        raise SystemExit("no member of %.1f-%.1f MB found before the stop offset" % (
            options.min_mb, options.max_mb))
    return saved, found


def compare_with_manifest(saved):
    """Extracted members against the committed hash list. Returns the differences.
    Entries of the list that were not extracted now (other files kept in the corpus,
    such as the legacy prefix below) are not differences."""
    if not HASH_MANIFEST.exists():
        return ["new corpus (no committed hash list yet)"]
    old = {}
    for line in HASH_MANIFEST.read_text(encoding="utf-8").splitlines():
        if line.strip() and not line.startswith("#"):
            digest, relative = line.split(None, 1)
            old[relative.strip()] = digest
    differences = []
    for relative in sorted(saved):
        if relative not in old:
            differences.append("NEW " + relative)
        elif old[relative] != saved[relative][1]:
            differences.append("DIFFERS " + relative)
    return differences


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--from-file", type=Path,
                        help="parse a local prefix/download of the tar instead of the remote archive")
    parser.add_argument("--stop-after-mb", type=int, default=500)
    parser.add_argument("--members", default=",".join(PINNED),
                        help="comma-separated member base names to extract")
    parser.add_argument("--by-size", action="store_true",
                        help="select by --min-mb/--max-mb/--want instead of the member list")
    parser.add_argument("--want", type=int, default=8)
    parser.add_argument("--min-mb", type=float, default=0.5)
    parser.add_argument("--max-mb", type=float, default=8.0)
    parser.add_argument("--skip", type=int, default=0, help="members in the size window to skip")
    parser.add_argument("--update-manifest", action="store_true",
                        help="rewrite the hash list when the corpus differs from it")
    parser.add_argument("--verify", action="store_true",
                        help="only compare the corpus on disk with the committed hash list")
    options = parser.parse_args()
    if options.verify:
        return corpus_manifest.verify_hashes(HASH_MANIFEST)

    source = LocalSource(options.from_file) if options.from_file else RemoteSource(URL)
    saved, found = extract(source, options)
    differences = compare_with_manifest(saved)
    for line in differences:
        print("CORPUS " + line)
    if differences and HASH_MANIFEST.exists() and not options.update_manifest:
        print("the extracted corpus differs from the committed hash list; "
              "re-run with --update-manifest to accept the new bytes", file=sys.stderr)
        return 1

    CORPUS.mkdir(parents=True, exist_ok=True)
    for relative, (data, _) in saved.items():
        (ROOT / relative).write_bytes(data)
    notes = ["corpus of the eggNOG gate; member names are pinned in verify/eggnog_stream.py",
             "source: " + (options.from_file.as_posix() if options.from_file else URL)]
    if differences or not HASH_MANIFEST.exists():
        corpus_manifest.write_hashes(CORPUS, HASH_MANIFEST, notes)

    alignments = sorted(p for p in CORPUS.glob("*.faa"))
    rows = ["# eggNOG5 Bacteria clusters: %d alignments x %d methods" % (len(alignments), len(METHODS)),
            "# COG0037.trimmed_alg.prefix.faa is the first 12,083,206 bytes of COG0037: the"
            " truncated member an earlier extractor saved when the download window cut it."]
    for path in alignments:
        relative = path.relative_to(ROOT).as_posix()
        rows.extend("-in %s %s -out $OUT" % (relative, method) for method in METHODS)
    content = "\n".join(rows) + "\n"
    if not CASE_MANIFEST.exists() or CASE_MANIFEST.read_text(encoding="utf-8") != content:
        CASE_MANIFEST.write_text(content, encoding="utf-8", newline="\n")
        print("WROTE " + CASE_MANIFEST.relative_to(ROOT).as_posix())
    print("CASES %d" % (len(rows) - 2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
