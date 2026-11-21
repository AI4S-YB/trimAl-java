#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#
"""Confirm whether the acquired source matches the pinned v1.5.1 Git tree."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parent.parent
tree = json.loads((root / "verify/reference-tag-tree.json").read_text())
entries = {entry["path"]: entry["sha"] for entry in tree["tree"] if entry["type"] == "blob"}
source_files = json.loads((root / "verify/reference-source-sha256.json").read_text())
matched, different = [], []
for name in source_files:
    path = root / "upstream/source" / name
    data = path.read_bytes()
    blob = hashlib.sha1(b"blob " + str(len(data)).encode("ascii") + b"\0" + data).hexdigest()
    upstream = entries.get("source/" + name)
    (matched if blob == upstream else different).append(name)
report = {"tag": "v1.5.1", "commit": "d637091abe33595775f40480970d1a18d87a7bcb",
          "matched": matched, "different": different,
          "identical_source": not different}
(root / "verify/reference-tag-match.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
if not different:
    reference = root / "verify/reference.json"
    metadata = json.loads(reference.read_text())
    metadata.update({"pinned_tag": report["tag"], "pinned_commit": report["commit"],
                     "note": "All 24 acquired C++ source/header/makefile blobs match the v1.5.1 tag."})
    reference.write_text(json.dumps(metadata, indent=2) + "\n")
