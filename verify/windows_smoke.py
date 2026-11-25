#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#
"""Compare the final JAR on Windows Java with the same Linux oracle and input paths."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parent.parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--java-exe", default=os.environ.get("TRIMAL_WINDOWS_JAVA"))
    parser.add_argument("--jar", help="JAR to test (default: java/target/trimal.jar)")
    options = parser.parse_args()
    executable = options.java_exe
    if executable is None:
        result = subprocess.run(["bash", "-c", "command -v java.exe"], capture_output=True, text=True, check=True)
        executable = result.stdout.strip()
    if not executable:
        raise SystemExit("set TRIMAL_WINDOWS_JAVA to the Windows java.exe path reachable from WSL")
    path = Path(options.jar) if options.jar else ROOT / "java/target/trimal.jar"
    jar_windows = subprocess.check_output(["wslpath", "-w", str(path)], text=True).strip()
    run_dir = ROOT / "verify/windows-main" / ("run-" + time.strftime("%Y%m%d-%H%M%S"))
    output = run_dir / "cases"
    output.mkdir(parents=True)
    commands = [
        [], ["--version"], ["-h"], ["--unknown"],
        ["-in", "verify/corpus/generated/protein.fasta", "-backtrans", "verify/corpus/generated/coding.fasta", "-fasta"],
    ]
    files = ["upstream/dataset/example.004.AA.fasta",
             "verify/corpus/generated/random_03.fasta",
             "verify/corpus/generated/random_07.fasta"]
    flags = [["-gappyout"], ["-strict"], ["-strictplus"], ["-automated1"],
             ["-gt", "0.8", "-cons", "65", "-w", "2"], ["-noallgaps"],
             ["-selectcols", "{", "0,2,4-7", "}"], ["-sident", "-soverlap"],
             ["-sgc", "-sgt", "-ssc", "-sst"], ["-maxidentity", "0.7"],
             ["-clustal"], ["-nexus"], ["-phylip3.2"], ["-fasta_m10", "-keepheader"]]
    for file in files:
        commands.extend([["-in", file] + flag for flag in flags])
    records = []
    for index, args in enumerate(commands, 1):
        original = subprocess.run([str(ROOT / "verify/bin/trimal")] + args, cwd=ROOT,
                                  capture_output=True, timeout=30)
        java = subprocess.run([executable, "-XX:-UsePerfData", "-cp", jar_windows, "trimal.TrimAl"] + args,
                              cwd=ROOT, capture_output=True, timeout=30)
        directory = output / ("case_%03d" % index)
        directory.mkdir(exist_ok=True)
        for prefix, result in (("orig", original), ("java", java)):
            (directory / (prefix + ".stdout")).write_bytes(result.stdout)
            (directory / (prefix + ".stderr")).write_bytes(result.stderr)
        passed = (original.returncode == java.returncode and original.stdout == java.stdout
                  and original.stderr == java.stderr)
        records.append({"case": index, "args": args, "passed": passed,
                        "orig_rc": original.returncode, "java_rc": java.returncode})
    summary = {"total": len(records), "passed": sum(r["passed"] for r in records),
               "failed": sum(not r["passed"] for r in records)}
    report = {"suite": "Windows main CLI vs Linux reference", "java": executable,
              "summary": summary, "cases": records}
    (run_dir / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(summary))
    print("REPORT " + str(run_dir / "result.json"))
    return int(summary["failed"] > 0)


if __name__ == "__main__":
    raise SystemExit(main())
