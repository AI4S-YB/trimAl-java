#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#
"""Print completed differential results and the first differing stream bytes."""
from pathlib import Path
import collections
import json
import sys

work = Path(sys.argv[1])
rows = [json.loads(path.read_text()) for path in sorted(work.glob("*/result.json"))]
print("finished=%d counts=%s" % (len(rows), dict(collections.Counter(r["status"] for r in rows))))
for row in rows:
    if row["status"] not in ("FAIL", "SKIP", "KNOWN_DIVERGENCE"):
        continue
    print(row["status"], row["case"], row["differences"], row["args"])
    if row["status"] == "SKIP":
        continue
    directory = work / row["case"]
    for stream in ("stdout", "stderr"):
        original = (directory / ("orig." + stream)).read_bytes()
        java = (directory / ("java." + stream)).read_bytes()
        if original == java:
            continue
        end = min(len(original), len(java))
        first = next((i for i in range(end) if original[i] != java[i]), end)
        print("  %s bytes=%d/%d first=%d orig=%r java=%r" %
              (stream, len(original), len(java), first,
               original[max(0, first - 60):first + 160],
               java[max(0, first - 60):first + 160]))
