#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Compare what the reference oracle itself printed in two environments.

    python3 verify/numeric/compare_env.py RUN_DIR_A RUN_DIR_B [--write OUT.json]

Each run directory comes from run_numeric.py (cases.txt, cpp.txt) or run_math.py
(cases.txt, cpp.txt) in one build environment: `cpp.txt` is the C++ oracle's own
output, the port is not involved. The comparison shows in which fields, and how
often, the reference's float bits depend on the environment they were produced in.
"""
import argparse
import collections
import json
from pathlib import Path
import sys


def rows(path):
    """{case id: {field: value}} of one oracle output."""
    result = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "|" not in line:
            continue
        identifier, fields = line.split("|", 1)
        result[identifier] = dict(item.split("=", 1) for item in fields.split("|") if "=" in item)
    return result


def kinds(cases):
    """{case id: kind letter} of the case list."""
    result = {}
    for line in cases.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) >= 2:
            result[parts[1]] = parts[0]
    return result


def compare(directory_a, directory_b):
    left, right = rows(directory_a / "cpp.txt"), rows(directory_b / "cpp.txt")
    kind_of = kinds(directory_a / "cases.txt")
    shared = sorted(set(left) & set(right))
    per_kind = collections.defaultdict(lambda: {"cases": 0, "differing_cases": 0, "fields": collections.Counter()})
    for identifier in shared:
        kind = kind_of.get(identifier, "?")
        entry = per_kind[kind]
        entry["cases"] += 1
        differing = [field for field in set(left[identifier]) | set(right[identifier])
                     if left[identifier].get(field) != right[identifier].get(field)]
        if differing:
            entry["differing_cases"] += 1
            for field in differing:
                entry["fields"][field] += 1
    examples = []
    for identifier in shared:
        differing = sorted(field for field in set(left[identifier]) | set(right[identifier])
                           if left[identifier].get(field) != right[identifier].get(field))
        if differing and len(examples) < 12:
            first = differing[0]
            examples.append({"id": identifier, "kind": kind_of.get(identifier), "field": first,
                             "a": left[identifier].get(first), "b": right[identifier].get(first)})
    return {"cases_compared": len(shared), "only_in_a": len(set(left) - set(right)),
            "only_in_b": len(set(right) - set(left)),
            "by_kind": {kind: {"cases": entry["cases"], "differing_cases": entry["differing_cases"],
                               "fields": dict(entry["fields"])}
                        for kind, entry in sorted(per_kind.items())},
            "examples": examples}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory_a", type=Path)
    parser.add_argument("directory_b", type=Path)
    parser.add_argument("--write", type=Path)
    options = parser.parse_args()
    result = compare(options.directory_a, options.directory_b)
    if options.write:
        options.write.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    print("cases compared: %d (only in A: %d, only in B: %d)" % (
        result["cases_compared"], result["only_in_a"], result["only_in_b"]))
    for kind, entry in result["by_kind"].items():
        print("  %-2s %6d cases, %6d differ  %s" % (
            kind, entry["cases"], entry["differing_cases"],
            ", ".join("%s=%d" % item for item in sorted(entry["fields"].items())) or "-"))
    for example in result["examples"]:
        print("  e.g. %s %s %s: %s vs %s" % (example["kind"], example["id"], example["field"],
                                            example["a"], example["b"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
