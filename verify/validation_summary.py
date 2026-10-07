#!/usr/bin/env python3
"""Persist acceptance counts from completed reports without rounding away divergences."""
import argparse
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent.parent


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--main-report", type=Path)
    options = parser.parse_args()
    if options.main_report:
        main_report = options.main_report
    else:
        reports = list((ROOT / "verify/results").glob("quick_manifest-*/report.json"))
        if not reports:
            raise SystemExit("no completed integrated main report")
        main_report = max(reports, key=lambda p: p.stat().st_mtime)
    main_data = load(main_report)
    companions = load(ROOT / "verify/companions/result.json")
    windows = load(ROOT / "verify/companions/result_windows.json")
    windows_main = load(ROOT / "verify/windows-main/result.json")
    numeric = load(ROOT / "verify/numeric/results.json")
    math = load(ROOT / "verify/numeric/math.results.json")
    audit_path = ROOT / "verify/alignment_audit/results_final-correct/summary.json"
    if not audit_path.exists():
        raise SystemExit("final corrected boundary/round-trip audit is not complete")
    audit = load(audit_path)

    def latest_report(pattern):
        reports = sorted((ROOT / "verify/results").glob(pattern), key=lambda p: p.stat().st_mtime)
        if not reports:
            raise SystemExit("missing differential report: " + pattern)
        return load(reports[-1])

    extended = latest_report("extended_manifest-*/report.json")
    legacy = latest_report("legacy_manifest-*/report.json")
    scale = latest_report("scale_manifest-*/report.json")
    jar = ROOT / "java/target/trimal.jar"
    tests = {"tests": 0, "failures": 0, "errors": 0, "skipped": 0}
    for path in (ROOT / "java/target/surefire-reports").glob("TEST-*.xml"):
        element = ET.parse(path).getroot()
        for key in tests:
            tests[key] += int(element.attrib.get(key, 0))
    summary = {
        "reference": load(ROOT / "verify/reference.json"),
        "jar": {"path": str(jar.relative_to(ROOT)), "bytes": jar.stat().st_size,
                "sha256": hashlib.sha256(jar.read_bytes()).hexdigest(), "release": 11},
        "unit_tests": tests,
        "main_cli": main_data["summary"],
        "extended_cli": extended["summary"],
        "legacy_cli": legacy["summary"],
        "scale_cli": scale["summary"],
        "companions": companions["summary"],
        "windows_companions": windows["summary"],
        "windows_main_cli": windows_main["summary"],
        "numeric_modules": {k: numeric[k] for k in ("counts", "matches", "total", "exact_matches", "exit_codes")},
        "float_math": {k: math[k] for k in ("counts", "matches", "total", "exact_matches", "exit_codes")},
        "alignment_audit": {k: audit[k] for k in ("count", "match", "known_native_ub",
                            "native_crash_safety_divergence", "actual_fail", "timeouts", "java_exceptions",
                            "roundtrip_requested", "roundtrip_executed")},
        "scale_timings": {row["case"]: {"args": row["args"][0] + " " + row["args"][1] + " " + row["args"][2],
                                        "orig_seconds": row["orig"]["seconds"],
                                        "java_seconds": row["java"]["seconds"]}
                          for row in scale["cases"]},
        "reports": {"main_cli": str(main_report.resolve().relative_to(ROOT)),
                    "companions": "verify/companions/result.json",
                    "windows": "verify/companions/result_windows.json",
                    "numeric": "verify/numeric/results.json",
                    "math": "verify/numeric/math.results.json",
                    "alignment": str(audit_path.relative_to(ROOT))},
    }
    summary["gated_failures"] = (tests["failures"] + tests["errors"]
                                 + main_data["summary"]["fail"]
                                 + extended["summary"]["fail"]
                                 + legacy["summary"]["fail"]
                                 + scale["summary"]["fail"]
                                 + companions["summary"]["failed"]
                                 + windows["summary"]["failed"]
                                 + windows_main["summary"]["failed"]
                                 + len(numeric.get("differences", []))
                                 + len(math.get("differences", []))
                                 + audit["actual_fail"] + audit["timeouts"])
    output = ROOT / "verify/validation-summary.json"
    output.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return int(summary["gated_failures"] > 0)


if __name__ == "__main__":
    raise SystemExit(main())
