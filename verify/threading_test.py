#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#
"""Thread-count determinism gate.

For conservation/HTML-heavy workloads, run the Java build with several
-Dtrimal.threads values (including 1 = plain serial loop) and require every
stdout/stderr/exit status/output file to be byte-identical to the reference
binary. Repeated default-thread runs additionally guard against scheduling
nondeterminism. Run on Linux/WSL: python3 verify/threading_test.py [--jar <build>]
"""
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "verify"))
OUT = ROOT / "verify/threading"
CONFIGS = ["1", "2", "7", "0"]  # 0 = JVM default (no property)
REPEATS_DEFAULT = 2
# the JAR under test; `--jar <path>` selects another build
jar = str(ROOT / "java/target/trimal.jar")

import toolchain  # noqa: E402


def run(command, timeout):
    started = time.monotonic()
    process = subprocess.Popen(command, cwd=ROOT, stdin=subprocess.DEVNULL,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               start_new_session=True)
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, 9)
        stdout, stderr = process.communicate()
        return {"rc": 124, "stdout": stdout, "stderr": stderr, "seconds": 0.0, "timeout": True}
    return {"rc": process.returncode, "stdout": stdout, "stderr": stderr,
            "seconds": round(time.monotonic() - started, 1), "timeout": False}


def snapshot(directory):
    return {p.relative_to(directory).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(directory.rglob("*")) if p.is_file()}


def main():
    global jar
    if "--jar" in sys.argv:
        jar = sys.argv[sys.argv.index("--jar") + 1]
    native_bin = str(ROOT / "verify/bin/trimal")
    work = OUT / ("run-" + time.strftime("%Y%m%d-%H%M%S"))
    work.mkdir(parents=True, exist_ok=True)

    # every workload touches a parallelized loop (conservation, identity, spurious,
    # HTML coloring) while staying inside a sane time budget
    workloads = [
        ("strict_mid", ["-in", "upstream/dataset/example.015.AA.bctoNOG.ENOG41099F3.fasta",
                        "-strict", "-out", "$OUT"], 600),
        ("html_mid", ["-in", "upstream/dataset/example.004.AA.fasta",
                      "-gappyout", "-htmlout", "$HTML", "-out", "$OUT"], 300),
        ("gaps_scale", ["-in", "upstream/dataset/example.014.AA.EggNOG.COG0591.fasta",
                        "-gappyout", "-out", "$OUT"], 600),
        ("spurious", ["-in", "upstream/dataset/example.004.AA.fasta",
                      "-resoverlap", "0.8", "-seqoverlap", "60", "-out", "$OUT"], 300),
    ]

    results = []
    for name, template, timeout in workloads:
        case_dir = work / name
        case_dir.mkdir(parents=True)
        runs = []
        configs = CONFIGS + ["0"] * REPEATS_DEFAULT
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            futures = {}
            for index, config in enumerate(configs):
                label = "threads_" + config if config != "0" else "default"
                if configs.index(config) != index:
                    label += "_repeat%d" % index  # default-thread repeats
                run_dir = case_dir / label
                native_dir = run_dir / "native"
                java_dir = run_dir / "java"
                native_dir.mkdir(parents=True, exist_ok=True)
                java_dir.mkdir(parents=True, exist_ok=True)
                native_args = [argument.replace("$OUT", str(native_dir / "o.aln"))
                               .replace("$HTML", str(native_dir / "o.html"))
                               for argument in template]
                java_args = [argument.replace("$OUT", str(java_dir / "o.aln"))
                             .replace("$HTML", str(java_dir / "o.html"))
                             for argument in template]
                native_result = run([native_bin] + native_args, timeout)
                command = [toolchain.java_tool("java"), "-XX:-UsePerfData"]
                if config != "0":
                    command += ["-Dtrimal.threads=" + config]
                command += ["-jar", jar] + java_args
                futures[index] = (config, label, native_dir, java_dir,
                                  pool.submit(run, command, timeout), native_result)
            runs = []
            for index, (config, label, native_dir, java_dir, future,
                        native_result) in sorted(futures.items()):
                result = future.result()
                java_files = snapshot(java_dir)
                native_files = snapshot(native_dir)
                identical = (result["rc"] == native_result["rc"]
                             and result["stdout"] == native_result["stdout"]
                             and result["stderr"] == native_result["stderr"]
                             and java_files == native_files
                             and not result["timeout"] and not native_result["timeout"])
                runs.append({"threads": config, "label": label,
                             "identical_to_native": identical,
                             "java_seconds": result["seconds"],
                             "native_seconds": native_result["seconds"],
                             "java_rc": result["rc"], "native_rc": native_result["rc"],
                             "java_files": java_files, "native_files": native_files})
        digests = {json.dumps(r["java_files"], sort_keys=True) for r in runs}
        consistent = len(digests) == 1
        passed = all(r["identical_to_native"] for r in runs) and consistent
        summary = {"workload": name, "passed": passed, "cross_thread_consistent": consistent,
                   "configs": [{"threads": r["threads"], "passed": r["identical_to_native"],
                                "java_seconds": r["java_seconds"]} for r in runs]}
        results.append(summary)
        print(json.dumps(summary), flush=True)
        if not passed:
            for r in runs:
                if not r["identical_to_native"]:
                    print("MISMATCH %s threads=%s rc=%s/%s" % (name, r["threads"],
                          r["java_rc"], r["native_rc"]), file=sys.stderr)

    report = {"suite": "thread-count determinism", "work": work.relative_to(ROOT).as_posix(),
              "results": results, "passed": all(r["passed"] for r in results)}
    (OUT / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    print("THREADING", "PASS" if report["passed"] else "FAIL")
    print("REPORT " + str(OUT / "results.json"))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
