#!/usr/bin/env python3
"""Stream-extract selected member alignments from a remote tar via HTTP Range.

Reads sequential chunks of the tarball, parses 512-byte headers, and saves members
whose size falls inside [min_mb, max_mb]. Members outside the size window are
skipped by advancing the remote offset without downloading their data. Stops after
`want` members or `stop_after_mb` of the archive.
"""
import argparse
import sys
import urllib.request
from pathlib import Path

URL = "http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2/2_trimmed_algs.tar"


def fetch_range(start, end, retries=5):
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(URL, headers={"Range": "bytes=%d-%d" % (start, end)})
            with urllib.request.urlopen(req, timeout=180) as resp:
                return resp.read()
        except Exception as e:  # transient proxy/9P/network errors
            last = e
            print("  fetch retry %d after error: %r" % (attempt + 1, e), file=sys.stderr)
    raise RuntimeError("range fetch failed: %r" % last)


def octal(raw):
    text = raw.rstrip(b"\0 ").decode("latin1")
    return int(text, 8) if text else 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stop-after-mb", type=int, default=500)
    parser.add_argument("--want", type=int, default=8)
    parser.add_argument("--min-mb", type=float, default=0.5)
    parser.add_argument("--max-mb", type=float, default=8.0)
    parser.add_argument("--out-dir", default="verify/corpus/eggnog")
    parser.add_argument("--skip", type=int, default=0, help="members to skip before extracting")
    options = parser.parse_args()

    out_dir = Path(options.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    chunk = 16 * 1024 * 1024
    offset = 0
    stop_at = options.stop_after_mb * 1024 * 1024
    lo = options.min_mb * 1024 * 1024
    hi = options.max_mb * 1024 * 1024
    extracted = []
    seen = 0

    while offset < stop_at and len(extracted) < options.want:
        buf = fetch_range(offset, min(offset + chunk, stop_at) - 1)
        if not buf:
            break
        pos = 0
        incomplete = False
        while pos + 512 <= len(buf):
            header = buf[pos:pos + 512]
            name = header[0:100].split(b"\0")[0].decode("latin1")
            if not name:
                break  # zero block = end of this chunk's parseable region
            size = octal(header[124:136])
            data_start = pos + 512
            padded_end = data_start + ((size + 511) // 512) * 512
            if padded_end > len(buf):
                offset += data_start  # next fetch resumes at this member's data
                incomplete = True
                break
            seen += 1
            if lo <= size <= hi and len(extracted) < options.want:
                if options.skip > 0 and seen <= options.skip:
                    print("  skip(early) %-38s %.1f MB" % (name, size / 1048576), flush=True)
                else:
                    target = out_dir / Path(name).name
                    if not target.exists():
                        target.write_bytes(buf[data_start:data_end])
                    extracted.append((name, size))
                    print("extracted %-45s %.1f MB" % (name, size / 1048576), flush=True)
            pos = padded_end
        if incomplete:
            continue
        offset += len(buf)

    print("SUMMARY consumed=%.0fMB extracted=%d seen=%d" %
          (offset / 1048576, len(extracted), seen), flush=True)
    for name, size in extracted:
        print("SAVED %s %d" % (name, size), flush=True)


if __name__ == "__main__":
    main()
