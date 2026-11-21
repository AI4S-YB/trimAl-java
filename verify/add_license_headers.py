#!/usr/bin/env python3
"""Add GPL license headers to all Java sources and compact notices to scripts.

Java files get the full block modeled on the upstream trimAl header style.
Scripts (shebang first) get a compact SPDX notice inserted after the shebang.
Idempotent: files that already carry the marker are skipped.
"""
import io
import os
from pathlib import Path

ROOT = Path("/mnt/j/ClaudeCodeDev/trimal")
MARKER = "SPDX-License-Identifier: GPL-3.0-or-later"

JAVA_HEADER = """/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.5.1 (Java 11 port): a tool for automated trimming of multiple
    sequence alignments, together with its readAl and statal companions.

    This file is part of the trimAl Java port, a derivative work of trimAl
    (https://github.com/inab/trimal).

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.

***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

"""

SCRIPT_HEADER = """#
# %s
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#
"""


def add_java_header(path):
    text = io.open(path, encoding="utf-8").read()
    if MARKER in text:
        return False
    io.open(path, "w", encoding="utf-8", newline="\n").write(JAVA_HEADER + text)
    return True


def add_script_notice(path):
    text = io.open(path, encoding="utf-8").read()
    if MARKER in text:
        return False
    lines = text.split("\n")
    notice = (SCRIPT_HEADER % MARKER).rstrip("\n").split("\n")
    if lines and lines[0].startswith("#!"):
        new = [lines[0]] + notice + lines[1:]
    else:
        new = notice + lines
    io.open(path, "w", encoding="utf-8", newline="\n").write("\n".join(new))
    return True


def main():
    changed, skipped = [], []
    for base in ["java/src/main/java", "java/src/test/java"]:
        for dirpath, _, files in os.walk(ROOT / base):
            for f in files:
                if f.endswith(".java"):
                    p = Path(dirpath) / f
                    (changed if add_java_header(p) else skipped).append(str(p.relative_to(ROOT)))
    script_files = [
        "verify/differential.py", "verify/generate_suite.py", "verify/eggnog_stream.py",
        "verify/pfam_extract.py", "verify/pffull_extract.py", "verify/inspect_results.py",
        "verify/validation_summary.py", "verify/threading_test.py", "verify/fuzz_suite.py",
        "verify/windows_matrix.py", "verify/windows_smoke.py", "verify/identify_reference.py",
        "verify/compare_trimal.sh", "verify/run_all.sh", "verify/pdl.sh",
        "java/trimal", "java/readal", "java/statal",
    ]
    for rel in script_files:
        p = ROOT / rel
        if p.exists():
            (changed if add_script_notice(p) else skipped).append(rel)
    print("headered %d files, skipped %d (already marked)" % (len(changed), len(skipped)))
    for rel in changed:
        print("  +", rel)


if __name__ == "__main__":
    main()
