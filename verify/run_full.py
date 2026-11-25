#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Run the verification gates in order and record one provenance-stamped result each.

    python3 verify/run_full.py                        # every gate, in the order below
    python3 verify/run_full.py --only numeric,math    # a subset, by name
    python3 verify/run_full.py --list                 # names and commands
    python3 verify/run_full.py --skip pffull          # everything except the named gates

Every gate is an existing runner. This driver only sequences the gates, keeps each
one's log and exit status, and extracts a small summary from the report the runner
names on its last `REPORT <path>` line. Before the first gate it records what is
being verified: the commit, whether the tree is clean, the JAR hash, the reference
binaries' hashes and versions, and the toolchain. The result of a run is
verify/work/run-<timestamp>/gates.json, which validation_summary.py turns into the
committed summary.

Gates that need data or a host this machine does not have are refused with the
reason, unless --skip-missing is given:

- pfam / pffull / eggnog / readal_big / statal_big need the downloaded corpora;
- the windows_* gates need a Windows java.exe reachable from WSL.

Timing claims come from this driver's full, sequential run: the heavy gates come
last, so their recorded seconds are not inflated by concurrent work.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
JAR = ROOT / "java" / "target" / "trimal.jar"
PYTHON = sys.executable or "python3"
# Per-case limits. The reference flushes its output line by line, which is slow on a
# network or 9P-mounted working directory and slower still when the machine is busy;
# a generous limit keeps a slow reference from failing a case without letting a hang
# through (a hang still fails, a few minutes later).
MATRIX_TIMEOUT = 300


def differential(manifest, timeout, jobs):
    return [PYTHON, "verify/differential.py", manifest, "--jar", "{jar}",
            "--timeout", str(timeout), "--jobs", str(jobs)]


def windows_matrix(manifest):
    return [PYTHON, "verify/windows_matrix.py", "--manifest", manifest, "--jar", "{jar}", "--jobs", "4"]


def companions():
    return [PYTHON, "verify/companions/verify_companions.py", "--java-jar", "{jar}"]


def corpus(name):
    return ROOT / "verify" / "corpus" / name


def fuzz(label):
    return [PYTHON, "verify/fuzz_suite.py", "--label", label, "--jar", "{jar}"]


def argfuzz(label):
    return [PYTHON, "verify/argfuzz_suite.py", "--label", label, "--jar", "{jar}"]


# name, command, what the gate needs on this machine. "{jar}" and "{stamp}" are filled
# in when the gate starts.
GATES = [
    ("generate", [PYTHON, "verify/generate_suite.py"], None),
    ("numeric", [PYTHON, "verify/numeric/run_numeric.py", "--label", "full"], None),
    ("math", [PYTHON, "verify/numeric/run_math.py"], None),
    ("ub", [PYTHON, "verify/numeric/run_ub.py"], None),
    ("numeric_junit", [PYTHON, "verify/numeric/run_junit.py"], "maven repository"),
    ("quick", differential("verify/quick_manifest.txt", MATRIX_TIMEOUT, 2), None),
    ("extended", differential("verify/extended_manifest.txt", MATRIX_TIMEOUT, 2), None),
    ("legacy", differential("verify/legacy_manifest.txt", MATRIX_TIMEOUT, 2), None),
    ("fuzz", fuzz("fuzz"), None),
    ("fuzz2", fuzz("fuzz2"), None),
    ("fuzz3", fuzz("fuzz3"), None),
    ("argfuzz", argfuzz("argfuzz"), None),
    ("argfuzz2", argfuzz("argfuzz2"), None),
    ("argfuzz3", argfuzz("argfuzz3"), None),
    ("companions", companions(), None),
    ("coverage", [PYTHON, "verify/coverage_map.py"], "upstream source"),
    ("audit", [PYTHON, "verify/alignment_audit/run_audit.py", "--round", "{stamp}"], None),
    ("probes", [PYTHON, "verify/alignment_audit/run_probes.py", "{stamp}"], None),
    ("threading", [PYTHON, "verify/threading_test.py", "--jar", "{jar}"], None),
    ("readal_big", differential("verify/readal_big_manifest.txt", 3600, 2), corpus("eggnog")),
    ("statal_big", differential("verify/statal_big_manifest.txt", 3600, 2), corpus("eggnog")),
    ("pfam", differential("verify/pfam_manifest.txt", 3600, 2), corpus("pfam")),
    ("eggnog", differential("verify/eggnog_manifest.txt", 3600, 1), corpus("eggnog")),
    ("windows_quick", windows_matrix("verify/quick_manifest.txt"), "windows java.exe"),
    ("windows_extended", windows_matrix("verify/extended_manifest.txt"), "windows java.exe"),
    ("windows_smoke", [PYTHON, "verify/windows_smoke.py", "--jar", "{jar}"], "windows java.exe"),
    ("windows_companions", [PYTHON, "verify/companions/verify_windows_companions.py",
                             "--java-jar", "{jar}"], "windows java.exe"),
    ("scale", differential("verify/scale_manifest.txt", 3600, 1), None),
    ("pffull", differential("verify/pffull_manifest.txt", 3600, 1), corpus("pffull")),
]
GATE_NAMES = [name for name, _, _ in GATES]


def sha256(path):
    sha = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            sha.update(block)
    return sha.hexdigest()


def load_average():
    """System load average when a gate starts and ends: another workload on the same
    machine inflates the recorded seconds, and the numbers should be able to show it."""
    try:
        return [round(value, 2) for value in os.getloadavg()]
    except OSError:  # not available on this platform
        return None


def capture(command):
    try:
        done = subprocess.run(command, cwd=ROOT, stdin=subprocess.DEVNULL,
                              capture_output=True, timeout=60)
        return {"rc": done.returncode,
                "stdout": done.stdout.decode("latin1").strip(),
                "stderr": done.stderr.decode("latin1").strip()}
    except (OSError, subprocess.TimeoutExpired) as error:
        return {"rc": None, "stdout": "", "stderr": repr(error)}


def git(*arguments):
    return subprocess.run(["git"] + list(arguments), cwd=ROOT, capture_output=True,
                          text=True).stdout.strip()


def has_windows_java():
    found = os.environ.get("TRIMAL_WINDOWS_JAVA")
    if found:
        return Path(found).exists()
    return bool(subprocess.run(["bash", "-c", "command -v java.exe"], capture_output=True,
                               text=True).stdout.strip())


def missing(requirement):
    """Why a gate cannot run here, or None when it can."""
    if requirement is None:
        return None
    if isinstance(requirement, Path):
        if not requirement.is_dir() or not any(requirement.iterdir()):
            return "corpus directory missing or empty: " + requirement.relative_to(ROOT).as_posix()
        return None
    if requirement == "upstream source":
        return None if (ROOT / "upstream/source/main.cpp").is_file() else "upstream/ source snapshot missing"
    if requirement == "windows java.exe":
        return None if has_windows_java() else "no Windows java.exe reachable (set TRIMAL_WINDOWS_JAVA)"
    if requirement == "maven repository":
        home = Path(os.environ.get("MAVEN_REPOSITORY") or Path.home() / ".m2" / "repository")
        return None if home.is_dir() else "no Maven repository for the JUnit jars"
    return "unknown requirement: " + requirement


def display(path):
    try:
        return Path(path).resolve().relative_to(ROOT).as_posix()
    except ValueError:  # outside the repository
        return Path(path).resolve().as_posix()


def provenance():
    dirty = [line for line in git("status", "--porcelain").splitlines() if line.strip()]
    binaries = {}
    for name in ("trimal", "readal", "statal"):
        path = ROOT / "verify" / "bin" / name
        binaries[name] = ({"sha256": sha256(path), "bytes": path.stat().st_size,
                           "--version": capture([str(path), "--version"])}
                          if path.is_file() else None)
    ldd = capture(["ldd", "--version"])
    return {
        "commit": git("rev-parse", "HEAD"),
        "dirty_files": dirty,
        "jar": ({"path": display(JAR), "sha256": sha256(JAR),
                 "bytes": JAR.stat().st_size} if JAR.is_file() else None),
        "reference_binaries": binaries,
        "environment": {
            "uname": capture(["uname", "-a"])["stdout"],
            "glibc": ldd["stdout"].splitlines()[0] if ldd["stdout"] else ldd["stderr"],
            "gpp": capture(["g++", "--version"])["stdout"].splitlines()[:1],
            "java": [line for line in capture(["java", "-version"])["stderr"].splitlines()][:3],
            "python": platform.python_version(),
            "cpu_count": os.cpu_count(),
            "loadavg_at_start": load_average(),
        },
    }


def summarize(report_path):
    """The compact result of a gate from the report its runner named."""
    if not report_path:
        return {}
    try:
        data = json.loads((ROOT / report_path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        return {"unreadable_report": repr(error)}
    if not isinstance(data, dict):
        return {}
    if isinstance(data.get("results"), list) and "passed" in data:
        # thread determinism: one entry per workload, each with its thread settings
        return {"passed": bool(data["passed"]), "workloads": len(data["results"]),
                "configs": sum(len(row.get("configs", [])) for row in data["results"])}
    if "summary" in data and isinstance(data["summary"], dict):
        summary = dict(data["summary"])
        for key in ("manifest", "work"):
            summary.pop(key, None)
        return summary
    keep = ("counts", "matches", "total", "exact_matches", "production_mismatches",
            "exit_codes", "passed", "count", "match", "known_native_ub",
            "native_crash_safety_divergence", "actual_fail", "timeouts", "java_exceptions",
            "roundtrip_requested", "roundtrip_executed", "differential_calls",
            "identity_alignments", "contracts_pass", "exit_code", "test_classes",
            "found", "started", "successful", "failed", "ported", "reachable", "gaps")
    return {key: data[key] for key in keep if key in data}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", help="comma-separated gate names")
    parser.add_argument("--skip", help="comma-separated gate names")
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--skip-missing", action="store_true",
                        help="leave out gates whose data or host this machine lacks")
    parser.add_argument("--jar", type=Path,
                        help="JAR to verify (default: java/target/trimal.jar); its hash is recorded")
    args = parser.parse_args()
    global JAR
    if args.jar:
        JAR = args.jar.resolve()
    if args.list:
        for name, command, requirement in GATES:
            print("%-22s %-16s %s" % (name, requirement or "-", " ".join(command)))
        return 0
    only = set(filter(None, args.only.split(","))) if args.only else None
    skip = set(filter(None, args.skip.split(","))) if args.skip else set()
    unknown = (only or set()) | skip
    unknown -= set(GATE_NAMES)
    if unknown:
        parser.error("unknown gate(s): " + ", ".join(sorted(unknown)))

    stamp = time.strftime("%Y%m%d-%H%M%S")
    work = ROOT / "verify" / "work" / ("run-" + stamp)
    # Concurrent runs started in the same second must not share a directory.
    serial = 0
    while True:
        try:
            work.mkdir(parents=True)
            break
        except FileExistsError:
            serial += 1
            work = ROOT / "verify" / "work" / ("run-%s-%d" % (stamp, serial))
    record = {"provenance": provenance(), "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
              "gates": [], "skipped": {}}
    print("RUN " + work.relative_to(ROOT).as_posix(), flush=True)
    print("COMMIT %s%s" % (record["provenance"]["commit"],
          "" if not record["provenance"]["dirty_files"] else
          " (dirty: %d file(s))" % len(record["provenance"]["dirty_files"])), flush=True)
    jar = record["provenance"]["jar"]
    print("JAR %s %s" % (jar["sha256"][:16], jar["path"]) if jar else "JAR missing", flush=True)
    if jar is None:
        print("no JAR: build it first with java/mvnw -f java/pom.xml package", file=sys.stderr)
        return 2
    fill = {"{stamp}": stamp, "{jar}": str(JAR)}

    failed = False
    for name, command, requirement in GATES:
        if (only and name not in only) or name in skip:
            continue
        reason = missing(requirement)
        if reason:
            if not args.skip_missing:
                print("MISSING %s: %s (use --skip-missing to leave it out)" % (name, reason),
                      file=sys.stderr)
                return 2
            record["skipped"][name] = reason
            print("SKIPPED %s (%s)" % (name, reason), flush=True)
            continue
        command = [fill.get(part, part) for part in command]
        log = work / (name + ".log")
        print("START %s" % name, flush=True)
        load_before = load_average()
        started = time.monotonic()
        with open(log, "wb") as handle:
            code = subprocess.call(command, cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT,
                                   stdin=subprocess.DEVNULL)
        seconds = round(time.monotonic() - started, 1)
        load_after = load_average()
        text = log.read_bytes().decode("latin1")
        reports = re.findall(r"^REPORT (.+)$", text, flags=re.MULTILINE)
        report = reports[-1].strip() if reports else None
        if report:
            try:
                report = Path(report).resolve().relative_to(ROOT).as_posix()
            except ValueError:
                pass
        entry = {"name": name, "command": command, "rc": code, "seconds": seconds,
                 "loadavg": [load_before, load_after],
                 "log": log.relative_to(ROOT).as_posix(), "report": report,
                 "summary": summarize(report)}
        record["gates"].append(entry)
        failed |= code != 0
        print("END %-16s rc=%-3s %8.1fs %s" % (
            name, code, seconds, json.dumps(entry["summary"], sort_keys=True)), flush=True)
        (work / "gates.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    record["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    (work / "gates.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print("GATES " + work.relative_to(ROOT).as_posix() + "/gates.json", flush=True)
    print("RESULT %s (%d gate(s), %d skipped)" % (
        "FAIL" if failed else "PASS", len(record["gates"]), len(record["skipped"])), flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
