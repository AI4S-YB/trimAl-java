#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Record what the reference executables print and write, for the in-process JUnit gate.

    python3 verify/record_golden.py --stats            # what would be recorded, and its size
    python3 verify/record_golden.py                    # write java/src/test/resources/golden/

Every case that the differential harness compared strictly (status PASS) in the
newest quick, extended and companions reports becomes one golden case: the reference's
exit status, stdout, stderr and produced files are stored, never the Java side's.
Known native divergences are left out: they are not matches. The JUnit test
GoldenCliTest then runs the port in-process on the same arguments and requires the
same bytes, on any operating system and JDK that can build the project, without the
reference binaries.

Layout under java/src/test/resources/golden/:

    root/           the input files at their repository-relative paths; the test runs
                    with this directory as working directory, so the arguments are the
                    very strings of the manifests (relative paths, `$OUT`, `$HTML`, `$IN`)
    expected.zip    per case: stdout, stderr and the produced files, as the
                    reference produced them
    cases.txt       the case list: id, tool, arguments, exit status, produced files

Cases that read an input larger than --max-input-bytes are left out: a handful of
multi-megabyte dataset files would dominate the repository for a few cases.
Only PASS cases are recorded, so a JUnit failure means the port changed, not that
the recording caught a moving target.
"""
import argparse
import json
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
GOLDEN = ROOT / "java/src/test/resources/golden"
RESULTS = ROOT / "verify/results"
DIFFERENTIAL = ["quick", "extended", "legacy"]
TOOLS = {"trimal": "trimal.TrimAl", "readal": "trimal.ReadAl", "statal": "trimal.StatAl"}
# options whose value is a file the command reads
INPUT_OPTIONS = ("-in", "-i", "-matrix", "-backtrans", "-forceselect", "-compareset")
FIXED_TIME = (1980, 1, 1, 0, 0, 0)  # zip entry timestamps, so re-recording is byte-stable


def newest(pattern):
    found = sorted(ROOT.glob(pattern), key=lambda p: p.stat().st_mtime)
    if not found:
        raise SystemExit("no report matching " + pattern)
    return found[-1]


def read(path):
    return path.read_bytes() if path.is_file() else None


def differential_cases(suite):
    """Golden cases of one differential manifest, from its newest report."""
    report_path = newest("verify/results/%s_manifest-*/report.json" % suite)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    work = report_path.parent
    tool = report["summary"].get("main", "trimal.TrimAl").split(".")[-1].lower()
    tool = {"trimal": "trimal", "readal": "readal", "statal": "statal"}[tool]
    cases = []
    for case in report["cases"]:
        if case["status"] != "PASS":
            continue
        directory = work / case["case"]
        expected = {"stdout": read(directory / "orig.stdout"), "stderr": read(directory / "orig.stderr")}
        for name in sorted(p.name for p in (directory / "orig").glob("*") if p.is_file()):
            expected[name] = read(directory / "orig" / name)
        cases.append({"id": "%s/%s" % (suite, case["case"]), "tool": tool, "args": case["args"],
                      "rc": case["orig"]["rc"], "expected": expected})
    return report_path, cases


def companions_cases():
    """Golden cases of the readal/statal suite, from its newest report."""
    report_path = newest("verify/companions/runs/*/result.json")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    cases = []
    for case in report["cases"]:
        if not case["passed"]:
            continue
        # Only the suite's own destination file becomes `$OUT`. Other -out values are
        # part of what the case tests (a directory, a path whose parent is missing).
        args, outputs = [], []
        tokens = list(case["args"])
        while tokens:
            token = tokens.pop(0)
            value = tokens.pop(0) if token == "-out" and tokens else None
            if value is not None and value.endswith("/output.aln"):
                outputs.append(value)
                args += ["-out", "$OUT"]
            elif value is not None:
                args += ["-out", value]
            else:
                args.append(token)
        directory = ROOT / case["artifacts"]
        expected = {"stdout": read(directory / "native.stdout"), "stderr": read(directory / "native.stderr")}
        if outputs:
            expected["o.aln"] = read(directory / "native.aln")
        cases.append({"id": "companions/%s_%s" % (case["tool"], case["name"].replace(".", "_")),
                      "tool": case["tool"], "args": args, "rc": case["native"]["rc"], "expected": expected})
    return report_path, cases


def referenced_inputs(args):
    """Repository-relative paths a case reads: option values and compare-list members."""
    paths, i = [], 0
    while i < len(args):
        token = args[i]
        if token in INPUT_OPTIONS and i + 1 < len(args):
            value = args[i + 1].replace("$IN/", "")
            paths.append(value)
            if token == "-compareset" and (ROOT / value).is_file():
                for line in (ROOT / value).read_bytes().split(b"\n"):
                    member = line.strip().replace(b"\r", b"").decode("latin1")
                    if member and (ROOT / member).is_file():
                        paths.append(member)
        i += 1
    return paths


def all_cases(selected):
    cases = []
    for suite in selected:
        if suite == "companions":
            report, found = companions_cases()
        else:
            report, found = differential_cases(suite)
        print("%-12s %s: %d strict matches" % (suite, report.parent.name, len(found)))
        cases += found
    return cases


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats", action="store_true", help="report what would be recorded")
    parser.add_argument("--suite", action="append",
                        help="differential manifest or 'companions' (default: quick, extended, companions)")
    parser.add_argument("--max-bytes", type=int, default=262144,
                        help="skip a case whose expected output is larger than this")
    parser.add_argument("--max-input-bytes", type=int, default=1000000,
                        help="skip a case that reads an input file larger than this")
    options = parser.parse_args()
    selected = options.suite or ["quick", "extended", "companions"]

    def largest_input(case):
        return max(((ROOT / relative).stat().st_size
                    for relative in referenced_inputs(case["args"]) if (ROOT / relative).is_file()),
                   default=0)

    cases, skipped, skipped_inputs, total = [], 0, 0, 0
    for case in all_cases(selected):
        if largest_input(case) > options.max_input_bytes:
            skipped_inputs += 1
            continue
        size = sum(len(data) for data in case["expected"].values() if data is not None)
        if size > options.max_bytes:
            skipped += 1
            continue
        total += size
        cases.append(case)

    inputs = {}
    for case in cases:
        for relative in referenced_inputs(case["args"]):
            if (ROOT / relative).is_file():
                inputs[relative] = (ROOT / relative).stat().st_size
    if "companions" in selected:
        for path in sorted((ROOT / "verify/companions/inputs").glob("*")):
            if path.is_file():
                inputs[path.relative_to(ROOT).as_posix()] = path.stat().st_size
    input_total = sum(inputs.values())
    print("cases %d (skipped %d over %d output bytes, %d over %d input bytes), "
          "expected bytes %d, input files %d (%d bytes)"
          % (len(cases), skipped, options.max_bytes, skipped_inputs, options.max_input_bytes,
             total, len(inputs), input_total))
    by_status = {}
    for case in cases:
        by_status[case["rc"]] = by_status.get(case["rc"], 0) + 1
    print("exit status of the reference among recorded cases:", sorted(by_status.items()))
    print("largest expected:", sorted(((sum(len(d) for d in c["expected"].values() if d), c["id"])
                                      for c in cases), reverse=True)[:5])
    absolute = [c["id"] for c in cases if any(b"/mnt/" in (d or b"") or b"C:\\\\" in (d or b"")
                                              for d in c["expected"].values())]
    print("cases whose expected output contains an absolute path:", absolute[:10], len(absolute))
    if options.stats:
        return 0

    # ---- write ------------------------------------------------------------------------
    root = GOLDEN / "root"
    if root.exists():
        for path in sorted(root.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            else:
                path.rmdir()
    for relative in sorted(inputs):
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((ROOT / relative).read_bytes())

    archive = GOLDEN / "expected.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zipped:
        for case in cases:
            for name, data in sorted(case["expected"].items()):
                if data is None:
                    continue
                info = zipfile.ZipInfo("%s/%s" % (case["id"], name), date_time=FIXED_TIME)
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = 0o644 << 16
                zipped.writestr(info, data)

    def escape(text):
        return text.replace("\\", "\\\\").replace("\n", "\\n").replace("\r", "\\r")

    reference = json.loads((ROOT / "verify/reference.json").read_text(encoding="utf-8"))
    lines = ["# Golden CLI cases of the trimAl Java port, recorded by verify/record_golden.py.",
             "# Reference: %s, %s %s (%s)." % (reference.get("identity", ""), reference.get("compiler", ""),
                                               reference.get("flags", ""), reference.get("pinned_commit", "")),
             "# Per case: id, tool and the reference's exit status; one line per argument; "
             "the files it produced.", ""]
    for case in cases:
        lines.append("case %s %s %d" % (case["id"], case["tool"], case["rc"]))
        lines += ["arg " + escape(argument) for argument in case["args"]]
        lines += ["file " + name for name in sorted(case["expected"])
                  if name not in ("stdout", "stderr") and case["expected"][name] is not None]
        lines.append("end")
    (GOLDEN / "cases.txt").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print("WROTE %s (%d cases, zip %d bytes, inputs %d bytes)"
          % (GOLDEN.relative_to(ROOT).as_posix(), len(cases), archive.stat().st_size, input_total))
    return 0


if __name__ == "__main__":
    sys.exit(main())
