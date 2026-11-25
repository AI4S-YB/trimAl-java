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
import math
from pathlib import Path
import struct
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


def ulps(a, b):
    """Distance in float32 ulps between two finite same-sign values, as bit patterns."""
    def position(bits):
        return bits if bits < 0x80000000 else 0x80000000 - (bits & 0x7FFFFFFF)
    return abs(position(a) - position(b))


def classify(left, right):
    """What kind of difference a field has: [(index, '1 ulp'), ...], one entry per
    differing element (a field may hold one value or a comma-separated vector)."""
    if left is None or right is None:
        return [(0, "missing")]
    left, right = left.split(","), right.split(",")
    if len(left) != len(right):
        return [(0, "length")]
    kinds = []
    for index, (a, b) in enumerate(zip(left, right)):
        if a == b:
            continue
        try:
            x, y = int(a), int(b)
        except ValueError:
            kinds.append((index, "text"))
            continue
        if not (0 <= x <= 0xFFFFFFFF and 0 <= y <= 0xFFFFFFFF):
            kind = "value"                 # an integer field, not a float bit pattern
        elif x ^ y == 0x80000000:
            kind = "sign bit only"         # e.g. the NaN returned for a negative input
        elif math.isnan(struct.unpack("<f", struct.pack("<I", x))[0]) or \
                math.isnan(struct.unpack("<f", struct.pack("<I", y))[0]):
            kind = "nan"
        else:
            distance = ulps(x, y)
            kind = ("%d ulp%s" % (distance, "" if distance == 1 else "s")
                    if distance <= 8 else "more than 8 ulps")
        kinds.append((index, kind))
    return kinds


def compare(directory_a, directory_b):
    left, right = rows(directory_a / "cpp.txt"), rows(directory_b / "cpp.txt")
    kind_of = kinds(directory_a / "cases.txt")
    shared = sorted(set(left) & set(right))
    per_kind = collections.defaultdict(lambda: {"cases": 0, "differing_cases": 0, "fields": collections.Counter()})
    difference_kinds = collections.defaultdict(collections.Counter)
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
                for _, item in classify(left[identifier].get(field), right[identifier].get(field)):
                    difference_kinds["%s %s" % (kind, field)][item] += 1
    # examples: prefer rounding differences over NaN-sign flips
    examples = []
    for identifier in shared:
        differing = sorted(field for field in set(left[identifier]) | set(right[identifier])
                           if left[identifier].get(field) != right[identifier].get(field))
        for field in differing:
            items = classify(left[identifier].get(field), right[identifier].get(field))
            if not any("ulp" in item for _, item in items):
                continue
            first, kind_of_difference = next((index, item) for index, item in items if "ulp" in item)
            a = left[identifier].get(field).split(",")
            b = right[identifier].get(field).split(",")
            if len(examples) < 12:
                examples.append({"id": identifier, "kind": kind_of.get(identifier), "field": field,
                                 "kind_of_difference": kind_of_difference,
                                 "a": a[first][:40], "b": b[first][:40]})
    return {"cases_compared": len(shared), "only_in_a": len(set(left) - set(right)),
            "only_in_b": len(set(right) - set(left)),
            "by_kind": {kind: {"cases": entry["cases"], "differing_cases": entry["differing_cases"],
                               "fields": dict(entry["fields"])}
                        for kind, entry in sorted(per_kind.items())},
            "difference_kinds": {field: dict(counts) for field, counts in sorted(difference_kinds.items())},
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
    for field, counts in result["difference_kinds"].items():
        print("     %-12s %s" % (field, ", ".join(
            "%d %s" % (count, item) for item, count in sorted(counts.items(), key=lambda item: -item[1]))))
    for example in result["examples"]:
        print("  e.g. %s %s %s (%s): %s vs %s" % (
            example["kind"], example["id"], example["field"], example.get("kind_of_difference", ""),
            example["a"], example["b"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
