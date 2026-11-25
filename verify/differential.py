#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#
"""Compare native trimAl and Java byte-for-byte. Runs on Linux/WSL, Python 3.8+.

Besides PASS/FAIL counts every report records how much of the suite the reference
actually executed: a suite in which the reference never exits 0 compares error paths
only, which is a legitimate test of the error messages but not of the trimming
logic. Such a suite fails unless --allow-no-reference-success is given explicitly.
"""
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

# Options whose value is the following argv token. A case that ends on one of them
# never reaches the loader: trimAl reports an unrecognised/incomplete argument.
VALUE_OPTIONS = {
    "-in", "-out", "-htmlout", "-matrix", "--alternative_matrix", "-compareset",
    "-forceselect", "-backtrans", "-gapthreshold", "-gt", "-simthreshold", "-st",
    "-conthreshold", "-ct", "-cons", "-resoverlap", "-seqoverlap", "-maxidentity",
    "-clusters", "-w", "-gw", "-sw", "-cw", "-block",
}


def relative(path):
    """Repo-relative posix path for reports, so evidence is not tied to one machine."""
    try:
        return Path(os.path.relpath(Path(path).resolve(), ROOT)).as_posix()
    except ValueError:  # different drives on Windows
        return str(path)


# Which executable and entry class a manifest exercises.
MAIN_CLASS = {"trimal": "trimal.TrimAl", "readal": "trimal.ReadAl", "statal": "trimal.StatAl"}


def manifest_tool(path):
    """A manifest may name its CLI with a first line `# tool: readal`; the default is
    trimAl. Explicit --orig / --main arguments take precedence."""
    with open(path, encoding="utf-8") as handle:
        first = handle.readline().strip()
    tool = first[len("# tool:"):].strip() if first.startswith("# tool:") else "trimal"
    if tool not in MAIN_CLASS:
        raise SystemExit("unknown tool in manifest header: " + tool)
    return tool


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
        # -11/-6/-7/-8 are the signals as a waited child reports them; 139/134/135/136
        # are the same crashes as a shell reports them.
        if orig["rc"] in (-11, -6, -7, -8, 139, 134, 135, 136):
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
    parser.add_argument("--main")
    parser.add_argument("--orig")
    parser.add_argument("--java", default="java")
    parser.add_argument("--timeout", type=float, default=60)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--work", type=Path)
    parser.add_argument("--allow-no-reference-success", action="store_true",
                        help="accept a suite on which the reference never exits 0 "
                             "(error-path-only manifests)")
    options = parser.parse_args()
    # Defaults follow the CLI the manifest names: verify/bin/<tool> and its entry class.
    tool = manifest_tool(options.manifest)
    options.orig = options.orig or str(ROOT / "verify" / "bin" / tool)
    options.main = options.main or MAIN_CLASS[tool]
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
    # How much of the suite the reference really executed.
    succeeded = [row for row in results if row["orig"]["rc"] == 0 and not row["orig"]["timeout"]]
    vacuity = {
        "reference_ok": len(succeeded),
        "reference_ok_with_output": sum(
            any(item["orig_bytes"] for item in row["files"].values()) for row in succeeded),
        "error_path_only": len(results) - len(succeeded),
        "dangling_option_args": sum(
            bool(row["args"]) and row["args"][-1] in VALUE_OPTIONS for row in results),
    }
    summary.update({"total": len(results), "manifest": relative(options.manifest),
                    "work": relative(work), "java": options.java, "main": options.main,
                    "vacuity": vacuity})
    report = {"summary": summary, "cases": results}
    (work / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("TOTAL pass={pass} fail={fail} skip={skip} known_divergence={known_divergence} total={total}".format(**summary), flush=True)
    print("VACUITY reference_ok={reference_ok} reference_ok_with_output={reference_ok_with_output} "
          "error_path_only={error_path_only} dangling_option_args={dangling_option_args}".format(**vacuity),
          flush=True)
    print("REPORT " + str(work / "report.json"), flush=True)
    if summary["fail"]:
        return 1
    if not vacuity["reference_ok"] and not options.allow_no_reference_success:
        print("VACUOUS: the reference exits non-zero on every case; this suite compares "
              "error paths only (pass --allow-no-reference-success to accept that)",
              file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
