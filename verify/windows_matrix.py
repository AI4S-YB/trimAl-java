#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#
"""Full cross-platform gate: run every quick-matrix case with Windows Java against
the Linux reference. Byte comparison of rc/stdout/stderr/files, no normalization.
Requires Windows java.exe reachable from WSL (TRIMAL_WINDOWS_JAVA or --java-exe).
Run inside WSL: python3 verify/windows_matrix.py [--jobs N]
"""
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "verify"))

NATIVE = str(ROOT / "verify/bin/trimal")


def wslpath_windows(p):
    return subprocess.check_output(["wslpath", "-w", str(p)], text=True).strip()


def run_side(command, out_dir, timeout):
    started = time.monotonic()
    process = subprocess.Popen(command, cwd=str(ROOT), stdin=subprocess.DEVNULL,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               start_new_session=True)
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, 9)
        stdout, stderr = process.communicate()
        return {"rc": 124, "stdout": stdout, "stderr": stderr, "timeout": True,
                "seconds": round(time.monotonic() - started, 1)}
    return {"rc": process.returncode, "stdout": stdout, "stderr": stderr,
            "timeout": False, "seconds": round(time.monotonic() - started, 1)}


def snapshot(directory):
    return {p.relative_to(directory).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(directory.rglob("*")) if p.is_file()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--java-exe", default=os.environ.get("TRIMAL_WINDOWS_JAVA"))
    parser.add_argument("--jar", default=str(ROOT / "java/target/trimal.jar"))
    parser.add_argument("--manifest", default=str(ROOT / "verify/quick_manifest.txt"))
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--filter", default="")
    options = parser.parse_args()

    executable = options.java_exe
    if executable is None:
        probe = subprocess.run(["bash", "-c", "command -v java.exe"], capture_output=True, text=True)
        executable = probe.stdout.strip()
    if not executable:
        raise SystemExit("Windows java.exe not found; set TRIMAL_WINDOWS_JAVA")

    jar_windows = wslpath_windows(Path(options.jar).resolve())
    manifest_path = Path(options.manifest)
    cases = []
    for raw in manifest_path.read_text(encoding="utf-8").splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("SKIP:"):
            continue
        cases.append(shlex.split(stripped))
    if options.filter:
        cases = [c for c in cases if options.filter in " ".join(c)]

    work = ROOT / "verify/windows-matrix" / time.strftime("%Y%m%d-%H%M%S")
    print("manifest=%s cases=%d" % (options.manifest, len(cases)), flush=True)
    work.mkdir(parents=True)
    results = []

    def one_case(index, arguments):
        name = "case_%04d" % index
        case_dir = work / name
        native_dir = case_dir / "native"
        windows_dir = case_dir / "windows"
        native_dir.mkdir(parents=True)
        windows_dir.mkdir(parents=True)
        native_args = [a.replace("$OUT", str(native_dir / "o.aln"))
                        .replace("$HTML", str(native_dir / "o.html"))
                        .replace("$IN", str(ROOT)) for a in arguments]
        windows_dir_win = "J:\\ClaudeCodeDev\\trimal\\verify\\windows-matrix\\" + work.name + "\\" + name + "\\windows"
        windows_args = [a.replace("$OUT", windows_dir_win + "\\o.aln")
                          .replace("$HTML", windows_dir_win + "\\o.html")
                          .replace("$IN", "J:\\ClaudeCodeDev\\trimal")
                        for a in arguments]
        native_result = run_side([NATIVE] + native_args, native_dir, options.timeout)
        command = [executable, "-XX:-UsePerfData", "-cp", jar_windows, "trimal.TrimAl"] + windows_args
        windows_result = run_side(command, windows_dir, options.timeout)
        (case_dir / "native.stdout").write_bytes(native_result["stdout"])
        (case_dir / "native.stderr").write_bytes(native_result["stderr"])
        (case_dir / "windows.stdout").write_bytes(windows_result["stdout"])
        (case_dir / "windows.stderr").write_bytes(windows_result["stderr"])
        java_files = snapshot(windows_dir)
        native_files = snapshot(native_dir)
        identical = (native_result["rc"] == windows_result["rc"]
                     and native_result["stdout"] == windows_result["stdout"]
                     and native_result["stderr"] == windows_result["stderr"]
                     and java_files == native_files
                     and not native_result["timeout"] and not windows_result["timeout"])
        row = {"case": name, "args": arguments, "status": "PASS" if identical else "FAIL",
               "native_rc": native_result["rc"], "windows_rc": windows_result["rc"],
               "native_seconds": native_result["seconds"], "windows_seconds": windows_result["seconds"]}
        if not identical:
            # Same evidence-based classification as differential.py: the only
            # difference is the out-of-bounds integer printCorrespondence prints for
            # a zero-column mapping (rc/stderr/files all match) -> known divergence.
            import re as _re
            if (native_result["rc"] == windows_result["rc"]
                    and native_result["stderr"] == windows_result["stderr"]
                    and java_files == native_files
                    and _re.fullmatch(rb"#ColumnsMap\t-?\d+\n", native_result["stdout"]) is not None
                    and windows_result["stdout"] == b"#ColumnsMap\t\n"):
                row["status"] = "KNOWN_DIVERGENCE"
                row["classification"] = "empty-map-oob-integer (native reads uninitialized heap)"
        (case_dir / "result.json").write_text(json.dumps(row, indent=2) + "\n", encoding="utf-8")
        return row

    with concurrent.futures.ThreadPoolExecutor(max_workers=options.jobs) as pool:
        futures = {pool.submit(one_case, i, c): i for i, c in enumerate(cases, 1)}
        done = 0
        for future in concurrent.futures.as_completed(futures):
            row = future.result()
            results.append(row)
            done += 1
            mark = row["status"]
            if mark != "PASS":
                print("%s %s native_rc=%s win_rc=%s" % (mark, row["case"],
                      row["native_rc"], row["windows_rc"]), flush=True)
            elif done % 50 == 0:
                print("... %d/%d" % (done, len(cases)), flush=True)

    summary = {"total": len(results),
               "passed": sum(r["status"] == "PASS" for r in results),
               "failed": sum(r["status"] == "FAIL" for r in results),
               "known_divergence": sum(r["status"] == "KNOWN_DIVERGENCE" for r in results)}
    report = {"suite": "Windows Java vs Linux reference, full quick matrix",
              "java_exe": executable, "jar": options.jar, "summary": summary,
              "cases": results}
    (work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(summary))
    print("REPORT " + str(work / "report.json"))
    return int(summary["failed"] > 0)


if __name__ == "__main__":
    sys.exit(main())
