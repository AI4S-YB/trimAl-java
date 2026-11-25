#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Confirm whether the acquired source matches the pinned v1.5.1 Git tree.

Every file listed in reference-source-sha256.json under upstream/source is hashed
as a Git blob and compared with the tag's tree (reference-tag-tree.json). The
result is written to reference-tag-match.json, and the exit status is non-zero when
any file differs, so CI can gate on it.
"""
import hashlib
import json
from pathlib import Path
import sys

root = Path(__file__).resolve().parent.parent
tree = json.loads((root / "verify/reference-tag-tree.json").read_text(encoding="utf-8"))
entries = {entry["path"]: entry["sha"] for entry in tree["tree"] if entry["type"] == "blob"}
source_files = json.loads((root / "verify/reference-source-sha256.json").read_text(encoding="utf-8"))
matched, different, absent = [], [], []
for name in source_files:
    path = root / "upstream/source" / name
    if not path.is_file():
        absent.append(name)
        continue
    data = path.read_bytes()
    blob = hashlib.sha1(b"blob " + str(len(data)).encode("ascii") + b"\0" + data).hexdigest()
    upstream = entries.get("source/" + name)
    (matched if blob == upstream else different).append(name)
report = {"tag": "v1.5.1", "commit": "d637091abe33595775f40480970d1a18d87a7bcb",
          "matched": matched, "different": different, "absent": absent,
          "identical_source": not different and not absent}
(root / "verify/reference-tag-match.json").write_text(
    json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
print(json.dumps({key: (len(value) if isinstance(value, list) else value)
                  for key, value in report.items()}, indent=2))
if different or absent:
    print("upstream/source does not match tag v1.5.1: different %s, absent %s" % (different, absent),
          file=sys.stderr)
    sys.exit(1)
reference = root / "verify/reference.json"
metadata = json.loads(reference.read_text(encoding="utf-8"))
metadata.update({"pinned_tag": report["tag"], "pinned_commit": report["commit"],
                 "note": "All %d acquired C++ source/header/makefile blobs match the v1.5.1 tag."
                         % len(matched)})
reference.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8", newline="\n")
