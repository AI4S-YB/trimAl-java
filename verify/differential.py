#!/usr/bin/env python3
"""Compare native trimAl and Java byte-for-byte. Runs on Linux/WSL, Python 3.8+."""
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent.parent


def run(command, timeout):
    started = time.monotonic()
    process = subprocess.Popen(command, cwd=ROOT, stdin=subprocess.DEVNULL,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               start_new_session=True)
    timed_out = False
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
    return {"rc": process.returncode, "stdout": stdout, "stderr": stderr,
            "timeout": timed_out, "seconds": round(time.monotonic() - started, 3)}


def read_manifest(path):
    cases = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        reason = None
        if stripped.startswith("SKIP:"):
            tag, stripped = stripped.split(None, 1)
            reason = tag[5:]
        # Deliberately disable shell comments so '#...' in input names is retained.
        cases.append({"args": shlex.split(stripped), "reason": reason})
    return cases


def one_case(index, case, work, options):
    name = "case_%04d" % index
    case_dir = work / name
    case_dir.mkdir(parents=True, exist_ok=False)
    row = {"case": name, "args": case["args"], "reason": case["reason"]}
    outputs = {}
    paths = {}
    for side in ("orig", "java"):
        directory = case_dir / side
        directory.mkdir()
        paths[side] = directory
        replacements = {"$IN": str(ROOT), "$OUT": str(directory / "o.aln"),
                        "$HTML": str(directory / "o.html")}
        arguments = []
        for argument in case["args"]:
            for key, value in replacements.items():
                argument = argument.replace(key, value)
            arguments.append(argument)
        if side == "orig":
            command = [options.orig] + arguments
        else:
            if options.classes:
                command = [options.java, "-XX:-UsePerfData", "-cp", options.classes, options.main] + arguments
            else:
                command = [options.java, "-XX:-UsePerfData", "-cp", options.jar, options.main] + arguments
        result = run(command, options.timeout)
        outputs[side] = result
        row[side] = {k: v for k, v in result.items() if k not in ("stdout", "stderr")}
        for stream in ("stdout", "stderr"):
            (case_dir / (side + "." + stream)).write_bytes(result[stream])
        (case_dir / (side + ".rc")).write_text(str(result["rc"]) + "\n", encoding="ascii")

    errors = []
    orig = outputs["orig"]
    java = outputs["java"]
    if orig["timeout"] or java["timeout"]:
        errors.append("timeout")
    if orig["rc"] != java["rc"]:
        errors.append("exit status")
    for stream in ("stdout", "stderr"):
        if orig[stream] != java[stream]:
            errors.append(stream)
    empty_map_ub = (orig["rc"] == java["rc"] and orig["stderr"] == java["stderr"]
                    and not orig["timeout"] and not java["timeout"]
                    and re.fullmatch(rb"#ColumnsMap\t-?\d+\n", orig["stdout"]) is not None
                    and java["stdout"] == b"#ColumnsMap\t\n")
    if empty_map_ub:
        # alignment.cpp:1807-1809 reads saveResidues[0] of a zero-length mapping when
        # a complementary/colnumbering run trims every column: the native integer is an
        # out-of-bounds read. rc, stderr and every produced file match; only that
        # integer differs. Classified, never counted as a byte-identical PASS.
        errors = [error for error in errors if error != "stdout"]
        if not errors:
            errors.append("empty-map-oob-integer")
    original_files = {p.relative_to(paths["orig"]): p for p in paths["orig"].rglob("*") if p.is_file()}
    java_files = {p.relative_to(paths["java"]): p for p in paths["java"].rglob("*") if p.is_file()}
    file_errors = []
    row["files"] = {}
    for relative in sorted(set(original_files) | set(java_files)):
        o = original_files.get(relative)
        j = java_files.get(relative)
        ob = o.read_bytes() if o else None
        jb = j.read_bytes() if j else None
        if ob != jb:
            file_errors.append(str(relative))
        row["files"][str(relative)] = {
            "orig_bytes": len(ob) if ob is not None else None,
            "java_bytes": len(jb) if jb is not None else None,
            "identical": ob == jb,
            "orig_sha256": hashlib.sha256(ob).hexdigest() if ob is not None else None,
            "java_sha256": hashlib.sha256(jb).hexdigest() if jb is not None else None,
        }
    errors.extend("file:" + value for value in file_errors)
    reason = case["reason"]
    row["status"] = "PASS" if not errors else "FAIL"
    # A known native crash is NOT a blanket permission to hide arbitrary output errors.
    if errors and reason == "orig-segv" and not orig["timeout"] and not java["timeout"]:
        if orig["rc"] in (-11, -6, 11, 139, 134):
            row["status"] = "KNOWN_DIVERGENCE"
            row["classification"] = "reference crashed; Java uses defined memory-safe behavior"
    elif errors and "empty-map-oob-integer" in errors:
        row["status"] = "KNOWN_DIVERGENCE"
        row["classification"] = ("only difference is the out-of-bounds integer printCorrespondence "
                                 "reads for a zero-column mapping; rc/stderr/files all match")
    elif errors and reason == "orig-ub-garbage-columnsmap":
        prefix = b"#ColumnsMap\t"
        if (orig["rc"] == java["rc"] and orig["stderr"] == java["stderr"]
                and not file_errors and not orig["timeout"] and not java["timeout"]
                and re.fullmatch(rb"#ColumnsMap\t-?\d+\n", orig["stdout"])
                and java["stdout"] == prefix + b"\n"):
            row["status"] = "SKIP"
    row["differences"] = errors
    (case_dir / "result.json").write_text(json.dumps(row, indent=2) + "\n", encoding="utf-8")
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--jar", default=str(ROOT / "java/target/trimal.jar"))
    parser.add_argument("--classes")
    parser.add_argument("--main", default="trimal.TrimAl")
    parser.add_argument("--orig", default=str(ROOT / "verify/bin/trimal"))
    parser.add_argument("--java", default="java")
    parser.add_argument("--timeout", type=float, default=60)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--work", type=Path)
    options = parser.parse_args()
    cases = read_manifest(options.manifest.resolve())
    if not cases:
        raise SystemExit("manifest contains no cases")
    work = (options.work or ROOT / "verify/results" /
            (options.manifest.stem + "-" + time.strftime("%Y%m%d-%H%M%S"))).resolve()
    work.mkdir(parents=True, exist_ok=False)
    if options.classes:
        options.classes = str(Path(options.classes).resolve())
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=options.jobs) as executor:
        pending = [executor.submit(one_case, i, case, work, options)
                   for i, case in enumerate(cases, 1)]
        for future in concurrent.futures.as_completed(pending):
            result = future.result()
            results.append(result)
            print("%s %s %s" % (result["status"], result["case"],
                                 ", ".join(result["differences"])), flush=True)
    results.sort(key=lambda item: item["case"])
    summary = {status.lower(): sum(row["status"] == status for row in results)
               for status in ("PASS", "FAIL", "SKIP", "KNOWN_DIVERGENCE")}
    summary.update({"total": len(results), "manifest": str(options.manifest.resolve()),
                    "work": str(work), "java": options.java, "main": options.main})
    report = {"summary": summary, "cases": results}
    (work / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("TOTAL pass={pass} fail={fail} skip={skip} known_divergence={known_divergence} total={total}".format(**summary), flush=True)
    print("REPORT " + str(work / "report.json"), flush=True)
    return 1 if summary["fail"] else 0


if __name__ == "__main__":
    sys.exit(main())
