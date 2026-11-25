#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Turn run_full.py records into the committed summary and the README tables.

    python3 verify/validation_summary.py                  # newest complete run
    python3 verify/validation_summary.py --run DIR [--run DIR ...]
    python3 verify/validation_summary.py --allow-partial  # accept a run that lacks gates
    python3 verify/validation_summary.py --render         # rewrite the README tables

A record is verify/work/run-<timestamp>/gates.json: one entry per gate with its exit
status, seconds, load, report and a compact summary. Several records may be merged
when they verified the same build - the same JAR and the same reference binaries -
even from different commits (a corpus changes, and only the gates that read it are
re-run); every gate then names the commit and the run it comes from, and a run made
in a dirty tree is refused, so the summary always describes one exact build.

A known native divergence is never counted as a match. `gated_failures` is the number
of gates with a non-zero exit status plus every failed case inside the other gates;
it must be zero for a summary to be written at all.

The README tables live between `<!-- verification-table:start -->` / `:end -->`,
`<!-- verification-provenance:... -->`, `<!-- scale-table:... -->` and
`<!-- mega-table:... -->` markers; `--render` replaces only what is between them, so
no number in those tables is typed by hand.
"""
import argparse
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "verify"))

import run_full  # noqa: E402  (GATE_NAMES)

WORK = ROOT / "verify" / "work"
OUTPUT = ROOT / "verify" / "validation-summary.json"
README = ROOT / "README.md"

# label, gates, how to read the summary
ROWS = [
    ("Main `trimal` CLI — quick matrix", ["quick"]),
    ("Main `trimal` CLI — extended matrix", ["extended"]),
    ("Main `trimal` CLI — legacy baseline", ["legacy"]),
    ("Main `trimal` CLI — scale (3,583-sequence real data + 4×60,000-column wide)", ["scale"]),
    ("`readal` / `statal` CLI", ["companions"]),
    ("Boundary/alignment audit (parse, writers, edge shapes, round trips) ‡", ["audit"]),
    ("Direct-source probes (matrix order, locale, 96 identity alignments)", ["probes"]),
    ("Loader fuzz (three labels × 180 mutations)", ["fuzz", "fuzz2", "fuzz3"]),
    ("CLI argument fuzz (two random labels and one accepted-shape label, 400 vectors each)",
     ["argfuzz", "argfuzz2", "argfuzz3"]),
    ("Downloaded real data: Pfam-A.seed top-10 families", ["pfam"]),
    ("Downloaded real data: 8 eggNOG5 Bacteria clusters and the truncated head of a ninth "
     "(plus its unaligned remainder)", ["eggnog"]),
    ("Downloaded real data: Pfam-A.full mega-families (up to 366K sequences)", ["pffull"]),
    ("Companions on the eggNOG clusters (readal conversions / statal stats)",
     ["readal_big", "statal_big"]),
    ("Thread determinism (4 workloads × 6 thread settings)", ["threading"]),
    ("Windows: full quick matrix (Windows Java vs Linux reference)", ["windows_quick"]),
    ("Windows: full extended matrix (Windows Java vs Linux reference)", ["windows_extended"]),
    ("Windows: CLI smoke + companions (Windows Java vs Linux reference)",
     ["windows_smoke", "windows_companions"]),
    ("Numeric module oracle (bit-level, incl. parallel paths)", ["numeric"]),
    ("Float math oracle (bit-level: the float overloads and the two double pow() sites)",
     ["math"]),
    ("Native undefined-behaviour observations (repeated native runs vs Java)", ["ub"]),
    ("Reference coverage check (reachable methods and CLI options)", ["coverage"]),
    ("JUnit compatibility tests outside Maven", ["numeric_junit"]),
]


def load_run(directory):
    path = Path(directory)
    if path.is_dir():
        path = path / "gates.json"
    if not path.is_file():
        raise SystemExit("no gates.json under " + str(directory))
    return json.loads(path.read_text(encoding="utf-8"))


def fingerprints(provenance):
    return {name: info["sha256"] for name, info in sorted(provenance["reference_binaries"].items())
            if info}


def merge(directories, allow_partial, allow_dirty):
    """One record out of several run_full.py runs. They must have verified the same
    build - the same JAR and the same reference binaries - and none may come from a
    dirty tree. They may come from different commits: after a change to a corpus,
    say, only the gates that read it are re-run. Every gate then names the commit
    it ran on, and the record lists them all."""
    directories = sorted(Path(directory).resolve() for directory in directories)
    records = [load_run(directory) for directory in directories]
    base = dict(records[0]["provenance"])
    for record in records:
        other = record["provenance"]
        if other["jar"]["sha256"] != base["jar"]["sha256"]:
            raise SystemExit("runs to merge must share one JAR: %s vs %s" % (
                base["jar"]["sha256"][:16], other["jar"]["sha256"][:16]))
        if fingerprints(other) != fingerprints(base):
            raise SystemExit("runs to merge must share the reference binaries")
        if other["dirty_files"] and not allow_dirty:
            raise SystemExit("a run was made in a dirty tree (%d file(s)); a summary describes "
                             "one exact build (pass --allow-dirty to override)" % len(other["dirty_files"]))
    base["commits"] = sorted({record["provenance"]["commit"] for record in records})
    gates, skipped = {}, {}
    for directory, record in zip(directories, records):
        for entry in record["gates"]:
            gates[entry["name"]] = dict(entry, commit=record["provenance"]["commit"],
                                        run=directory.relative_to(ROOT).as_posix())
        skipped.update(record.get("skipped", {}))
    absent = [name for name in run_full.GATE_NAMES if name not in gates and name not in skipped]
    if absent and not allow_partial:
        raise SystemExit("the run(s) do not cover: " + ", ".join(absent))
    runs = [directory.relative_to(ROOT).as_posix() for directory in directories]
    return base, gates, skipped, absent, runs


def provenance_lines(provenance):
    """The paragraph under the table: what exactly was verified, and where."""
    jar, environment = provenance["jar"], provenance["environment"]
    binaries = ", ".join("%s `%s`" % (name, info["sha256"][:12])
                         for name, info in sorted(provenance["reference_binaries"].items()) if info)
    commits = provenance.get("commits") or [provenance["commit"]]
    fields = environment["uname"].split()
    system = fields[0] if fields else ""
    kernel = fields[2] if len(fields) > 2 else ""
    compiler = (environment.get("gpp") or [""])[0].split(",")[0]
    load = (environment.get("loadavg_at_start") or [None])[0]
    lines = [
        "**Provenance.** Recorded by `verify/run_full.py`. The full record is "
        "[verify/validation-summary.json](verify/validation-summary.json), where every gate "
        "names the commit and the run it comes from.",
        "",
        "- Build: commit%s %s; JAR SHA-256 `%s` (%s bytes)." % (
            "s" if len(commits) > 1 else "", ", ".join("`%s`" % commit[:10] for commit in commits),
            jar["sha256"][:16], num(jar["bytes"])),
        "- Reference: %s, built from the pinned source with %s." % (binaries, compiler or "g++"),
        "- Environment: %s %s, %s, %s." % (system, kernel, environment.get("glibc", ""),
                                           (environment.get("java") or [""])[0]),
    ]
    if load is not None:
        lines.append("- Load average of the machine before the first gate: %s on %s cores. The "
                     "machine was not idle, so the seconds in the record are those of a machine "
                     "under that load; each speed-up compares the two implementations under the "
                     "same load." % (load, environment.get("cpu_count")))
    return lines


def tally(name, summary):
    """(matched, compared, known differences, failed) of one gate."""
    if not summary:
        return None
    if name in ("numeric", "math"):
        return (summary.get("exact_matches"), summary.get("total"), 0,
                summary.get("production_mismatches"))
    if name == "audit":
        return (summary.get("match"), summary.get("count"),
                (summary.get("known_native_ub") or 0) + (summary.get("native_crash_safety_divergence") or 0),
                (summary.get("actual_fail") or 0) + (summary.get("timeouts") or 0)
                + (summary.get("java_exceptions") or 0))
    if name == "probes":
        return (summary.get("match"), summary.get("differential_calls"), 0,
                0 if summary.get("contracts_pass") else 1)
    if name == "ub":
        counts = summary.get("counts") or {}
        return ((counts.get("aa_exact_java_matches") or 0) + (counts.get("cut_run_exact_java_matches") or 0),
                (counts.get("aa_original_invocations") or 0) + (counts.get("cut_original_invocations") or 0),
                counts.get("known_ub_output_differences") or 0, 0)
    if name == "coverage":
        return (summary.get("ported"), summary.get("reachable"), 0, summary.get("gaps"))
    if name == "numeric_junit":
        return (summary.get("successful"), summary.get("found"), 0,
                (summary.get("failed") or 0) + (1 if summary.get("exit_code") else 0))
    if name == "threading":
        passed = bool(summary.get("passed"))
        return (summary.get("configs", 0) if passed else 0, summary.get("configs", 0), 0,
                0 if passed else 1)
    if "pass" in summary:      # differential.py
        return (summary.get("pass"), summary.get("total"),
                (summary.get("known_divergence") or 0) + (summary.get("skip") or 0),
                summary.get("fail") or 0)
    if "passed" in summary:    # companions, windows_matrix, windows_smoke
        return (summary.get("passed"), summary.get("total"), summary.get("known_divergence") or 0,
                summary.get("failed") or 0)
    return None


def acceptance(gates, names):
    """(accepted, rejected) counts the reference produced for a row's gates, when
    their reports carry the harness's vacuity block."""
    accepted = rejected = 0
    seen = False
    for name in names:
        vacuity = (gates.get(name) or {}).get("summary", {}).get("vacuity")
        if vacuity:
            seen = True
            accepted += vacuity["reference_ok"]
            rejected += vacuity["error_path_only"]
    return (accepted, rejected) if seen else None


def junit_counts():
    """Surefire totals, and the result of the golden CLI class on its own (one of its
    tests checks the stream redirection seam, the rest are recorded reference cases)."""
    tests = {"tests": 0, "failures": 0, "errors": 0, "skipped": 0, "golden": None}
    for path in sorted((ROOT / "java/target/surefire-reports").glob("TEST-*.xml")):
        element = ET.parse(path).getroot()
        counts = {key: int(element.attrib.get(key, 0)) for key in ("tests", "failures", "errors", "skipped")}
        for key in counts:
            tests[key] += counts[key]
        if element.attrib.get("name") == "trimal.GoldenCliTest":
            tests["golden"] = {"cases": counts["tests"] - 1,
                               "failed": counts["failures"] + counts["errors"]}
    return tests


def num(value):
    return "{:,}".format(value)


def verification_table(gates):
    lines = ["| Suite | Strict byte/bit matches | Separately recorded differences | Reference accepted / rejected |",
             "|---|---:|---|---|"]
    for label, names in ROWS:
        tallies = [tally(name, (gates.get(name) or {}).get("summary")) for name in names]
        if any(row is None for row in tallies):
            lines.append("| %s | not run | — | — |" % label)
            continue
        matched, compared, known, failed = (sum(row[index] or 0 for row in tallies) for index in range(4))
        if names == ["threading"]:
            matches = "all byte-identical" if not failed and matched else "MISMATCH"
        else:
            matches = "%s / %s" % (num(matched), num(compared))
        if failed:
            difference = "**%s FAILED**" % num(failed) + (
                ", %s native UB/crash path%s" % (num(known), "" if known == 1 else "s") if known else "")
        elif known:
            difference = "%s native UB/crash path%s" % (num(known), "" if known == 1 else "s")
        else:
            difference = "None"
        counts = acceptance(gates, names)
        lines.append("| %s | %s | %s | %s |" % (
            label, matches, difference, "%s / %s" % counts if counts else "—"))
    tests = junit_counts()
    golden = tests.get("golden")
    if golden:
        lines.append("| Golden CLI gate: recorded reference output, port run in-process in `mvn test` | %s / %s | %s | — |" % (
            num(golden["cases"] - golden["failed"]), num(golden["cases"]),
            "None" if not golden["failed"] else "**%s FAILED**" % num(golden["failed"])))
    if tests["tests"]:
        lines.append("| Integrated JUnit (`mvn package`), the golden gate included | %s / %s | %s | — |" % (
            num(tests["tests"] - tests["failures"] - tests["errors"]), num(tests["tests"]),
            "None" if not (tests["failures"] or tests["errors"]) else "**FAILED**"))
    else:
        lines.append("| Integrated JUnit (`mvn package`) | not run | — | — |")
    return lines


def timings(gate):
    """[(label, reference seconds, Java seconds)] from a differential report."""
    report = (gate or {}).get("report")
    if not report:
        return []
    data = json.loads((ROOT / report).read_text(encoding="utf-8"))
    rows = []
    for case in data["cases"]:
        args = case["args"]
        source = Path(args[1]).name
        for token in args[2:]:
            if token.startswith("-") and token not in ("-out", "-htmlout"):
                source += " " + token
        rows.append((source, case["orig"]["seconds"], case["java"]["seconds"]))
    return rows


def timing_table(rows, title):
    lines = ["| %s | Reference | Java | Speed-up |" % title, "|---|---:|---:|---:|"]
    reference_total = java_total = 0.0
    for label, reference, java in rows:
        reference_total += reference
        java_total += java
        lines.append("| `%s` | %.1f s | %.1f s | %.0f× |" % (label, reference, java, reference / max(java, 0.001)))
    lines.append("| **total** | **%.0f s** | **%.0f s** | **%.0f×** |" % (
        reference_total, java_total, reference_total / max(java_total, 0.001)))
    return lines


def mega_table(rows):
    """Every method of the largest family, then the total over all cases."""
    lines = ["| Family — method | Reference | Java | Speed-up |", "|---|---:|---:|---:|"]
    largest = rows[0][0].split(" ")[0] if rows else ""
    for label, reference, java in rows:
        if label.startswith(largest + " "):
            lines.append("| `%s` | %.0f s | %.0f s | %.0f× |" % (label, reference, java, reference / max(java, 0.001)))
    reference_total = sum(row[1] for row in rows)
    java_total = sum(row[2] for row in rows)
    lines.append("| all %d cases, up to 366K sequences | %.0f min | %.0f min | %.0f× |" % (
        len(rows), reference_total / 60, java_total / 60, reference_total / max(java_total, 0.001)))
    return lines


def replace_block(text, name, lines):
    start, end = "<!-- %s:start -->" % name, "<!-- %s:end -->" % name
    if start not in text or end not in text:
        raise SystemExit("README.md lacks the %s markers" % name)
    head, rest = text.split(start, 1)
    _, tail = rest.split(end, 1)
    return head + start + "\n" + "\n".join(lines) + "\n" + end + tail


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="append", type=Path, help="a run_full.py run directory")
    parser.add_argument("--allow-partial", action="store_true",
                        help="accept a run that lacks gates (never one with a failure)")
    parser.add_argument("--allow-dirty", action="store_true")
    parser.add_argument("--allow-failures", action="store_true",
                        help="write the summary even with gated failures (diagnosis only)")
    parser.add_argument("--render", action="store_true", help="rewrite the README tables")
    parser.add_argument("--no-write", action="store_true")
    options = parser.parse_args()
    # default: the newest record; several records must be named (and share one build)
    directories = options.run or [sorted(WORK.glob("run-*"))[-1]]
    if not options.run and not list(WORK.glob("run-*")):
        raise SystemExit("no run_full.py record under " + str(WORK))
    provenance, gates, skipped, absent, runs = merge(directories, options.allow_partial, options.allow_dirty)

    failed_gates = [name for name, gate in gates.items() if gate["rc"]]
    tallies = {name: tally(name, gate["summary"]) for name, gate in gates.items()}
    failed_cases = sum((row[3] or 0) for row in tallies.values() if row)
    tests = junit_counts()
    summary = {
        "provenance": provenance,
        "runs": runs,
        "gates": {name: {key: gate.get(key) for key in (
                    "command", "rc", "seconds", "loadavg", "commit", "run", "log", "report", "summary")}
                  for name, gate in sorted(gates.items())},
        "skipped_gates": skipped,
        "gates_not_run": absent,
        "unit_tests": tests,
        "scale_timings": timings(gates.get("scale")),
        "pffull_timings": timings(gates.get("pffull")),
        "failed_gates": failed_gates,
        "failed_cases": failed_cases,
        "gated_failures": len(failed_gates) + failed_cases + tests["failures"] + tests["errors"],
    }
    print(json.dumps({key: summary[key] for key in
                      ("runs", "failed_gates", "failed_cases", "gated_failures", "gates_not_run")}, indent=2))
    if summary["gated_failures"] and not options.allow_failures:
        # Nothing is written and the README is not rendered: a record with a failure
        # must not become the committed summary, however incomplete the run was allowed to be.
        print("refusing to write a summary with %d gated failure(s)" % summary["gated_failures"],
              file=sys.stderr)
        return 1
    if not options.no_write:
        OUTPUT.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    if options.render:
        text = README.read_text(encoding="utf-8")
        text = replace_block(text, "verification-table", verification_table(gates))
        text = replace_block(text, "verification-provenance", provenance_lines(provenance))
        if "scale" in gates:
            text = replace_block(text, "scale-table", timing_table(summary["scale_timings"], "Job"))
        if "pffull" in gates:
            text = replace_block(text, "mega-table", mega_table(summary["pffull_timings"]))
        README.write_text(text, encoding="utf-8", newline="\n")
        print("README tables rendered")
    print("SUMMARY " + OUTPUT.relative_to(ROOT).as_posix())
    return 1 if summary["gated_failures"] else 0


if __name__ == "__main__":
    sys.exit(main())
