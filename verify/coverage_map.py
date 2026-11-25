#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Machine check that the port covers every CLI option and every method of the reference.

    python3 verify/coverage_map.py              # write verify/coverage-map.json, fail on gaps
    python3 verify/coverage_map.py --print      # also print the findings

Two things are compared against upstream/source (tag v1.5.1):

1. Options. Every option string literal of main.cpp / readAl.cpp / statAl.cpp must
   occur in the Java entry class of the same tool, and the Java class must not accept
   an option the reference does not have.
2. Methods. Every method the reference defines *and calls* must have a method of the
   same name in the Java class that ports it. Two things are not gaps in the port and
   are reported separately:
   - a reference method that is defined but unreachable (nothing outside the
     unreachable set ever calls it): porting it could not change any output;
   - a helper that the port has under another name or inlined. Those are listed in
     PORTED_AS below with the Java symbol that took over, and the symbol is checked
     to exist, so the list cannot rot silently.
   Constructors, destructors and operator= have no method name to compare.

The method comparison is by name, so it proves nothing about behaviour; the
differential suites do that. It proves that nothing was forgotten.
"""
import argparse
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "upstream" / "source"
JAVA = ROOT / "java" / "src" / "main" / "java" / "trimal"
OUTPUT = ROOT / "verify" / "coverage-map.json"

# reference unit: (header, defining sources, Java files that port it)
UNITS = {
    "alignment": (["alignment.h"], ["alignment.cpp", "autAlignment.cpp", "rwAlignment.cpp"],
                  ["Alignment.java"]),
    "sequencesMatrix": (["sequencesMatrix.h"], ["sequencesMatrix.cpp"], ["stat/SequencesMatrix.java"]),
    "similarityMatrix": (["similarityMatrix.h"], ["similarityMatrix.cpp"], ["stat/SimilarityMatrix.java"]),
    "statisticsGaps": (["statisticsGaps.h"], ["statisticsGaps.cpp"], ["stat/StatisticsGaps.java"]),
    "statisticsConservation": (["statisticsConservation.h"], ["statisticsConservation.cpp"],
                               ["stat/StatisticsConservation.java"]),
    "statisticsFiles": (["statisticsFiles.h"], ["statisticsFiles.cpp"], ["StatAl.java"]),
    "compareFiles": (["compareFiles.h"], ["compareFiles.cpp"], ["CompareFiles.java"]),
    "utils": (["utils.h"], ["utils.cpp"], ["util/Utils.java", "util/CppFormat.java", "util/CppOut.java",
                                                   "util/Values.java"]),
}
# reference method -> (Java file, Java method) that implements it under another name
PORTED_AS = {
    "checkFile": ("Alignment.java", "formatInputAlignment"),      # ifstream open failure -> format 0
    "readLine": ("Alignment.java", "readLine"),                   # Cursor.readLine
    "readNumbers_StartEnd": ("util/Utils.java", "readNumbersStartEnd"),
    "swap": ("util/Utils.java", "quicksort"),                     # swaps are inlined
}
# reference units whose methods are all static and are called as `Class::name(`
STATIC_UNITS = {"utils", "compareFiles"}
TOOLS = {"trimal": ("main.cpp", "TrimAl.java"), "readal": ("readAl.cpp", "ReadAl.java"),
         "statal": ("statAl.cpp", "StatAl.java")}
CONTROL = {"for", "if", "while", "switch", "return", "else", "catch", "try", "do", "new",
           "throw", "case", "break", "continue", "assert", "synchronized"}
DECLARATION = re.compile(r"^[ ]{2,4}(?:(?:virtual|static|inline)\s+)?"
                         r"(?:[A-Za-z_][\w:<>]*[ &*]+)?(~?\w+)\s*\(")
# an out-of-class definition starts the line: `void alignment::clean(...) {`
DEFINITION_LINE = re.compile(r"^\s{0,2}[\w\s:&*<>]*\b\w+::(~?\w+)\s*\(")
UNQUALIFIED_CALL = re.compile(r"(?<![\w:])(~?\w+)\s*\(")
QUALIFIED_CALL = re.compile(r"(\w+)::(~?\w+)\s*\(")
JAVA_METHOD = re.compile(r"^\s*(?:(?:public|private|protected|static|final|synchronized|abstract)\s+)*"
                         r"(?:<[^>]+>\s+)?([A-Za-z_][\w.]*(?:<[^>]*>)?(?:\[\])*)\s+([A-Za-z_]\w*)\s*\(")
OPTION = re.compile(r'"(-{1,2}[A-Za-z_]\w*)"')


def text(path):
    return path.read_text(encoding="latin1")


def require(path):
    if not path.is_file():
        raise SystemExit("missing reference source: %s (clone upstream tag v1.5.1 first)" % path)
    return text(path)


def cpp_methods(header, sources, class_name):
    """(declared, defined) method names of one reference class."""
    declared, defined = set(), set()
    for line in text(header).splitlines():
        stripped = line.lstrip()
        if stripped.startswith(("*", "//", "/", "#")) \
                or not 2 <= len(line) - len(stripped) <= 4:   # class members only
            continue
        match = DECLARATION.match(line)
        if not match:
            continue
        name = match.group(1)
        if name in (class_name, "~" + class_name):
            continue
        declared.add(name)
    bodies = "\n".join(text(source) for source in sources)
    for name in declared:
        # a definition is `Type Class::name(` of this very class: `utils::name(` inside
        # an alignment method is a call to another class, not a definition
        if re.search(r"(?<![\w:])" + re.escape(class_name) + "::" + re.escape(name) + r"\s*\(", bodies):
            defined.add(name)
    return declared, defined


def call_graph():
    """[(enclosing reference method or None, unqualified names called, qualified names
    called)] for every line of the reference sources that is not a definition header.
    A top-level function starts at an unindented line that has a `(`: either an
    out-of-class definition `Type Class::method(` (its lines belong to `method`) or
    free code such as main() and the tools' helpers (owner None, where reachability
    starts)."""
    graph, owner = [], None
    for source in sorted(SOURCE.glob("*.cpp")):
        owner = None
        for line in text(source).splitlines():
            if line and not line[0].isspace() and not line.startswith(("//", "/*", "*", "#")) \
                    and "(" in line:
                definition = DEFINITION_LINE.match(line)
                owner = definition.group(1) if definition else None
                if definition:
                    continue
            graph.append((owner, set(UNQUALIFIED_CALL.findall(line)),
                          {qualified for _, qualified in QUALIFIED_CALL.findall(line)}))
    return graph


def reachable(methods, static_methods, graph):
    """Methods that entry-point code can reach. A member method is called by a bare
    `name(`; a static one also by `Class::name(`. Reachability is a fixpoint from the
    ownerless lines, so a method called only from unreachable code is unreachable too."""
    live, changed = set(), True
    while changed:
        changed = False
        for owner, bare, qualified in graph:
            if owner is not None and owner not in live:
                continue
            for name in methods:
                if name in live:
                    continue
                if name in bare or (name in static_methods and name in qualified):
                    live.add(name)
                    changed = True
    return live


def java_methods(files):
    names = set()
    for relative in files:
        for line in text(JAVA / relative).splitlines():
            stripped = line.lstrip()
            if stripped.startswith(("*", "//", "/", "#")):
                continue
            match = JAVA_METHOD.match(line)
            if match and match.group(1) not in CONTROL and match.group(2) not in CONTROL \
                    and "=" not in line.split(match.group(2) + "(", 1)[0]:
                names.add(match.group(2))
    return names


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--print", action="store_true")
    parser.add_argument("--no-write", action="store_true")
    options = parser.parse_args()

    result = {"options": {}, "methods": {}, "ported_as": {}, "gaps": 0}
    for tool, (cpp_entry, java_entry) in TOOLS.items():
        cpp = set(OPTION.findall(require(SOURCE / cpp_entry)))
        java = set(OPTION.findall(text(JAVA / java_entry)))
        entry = {"reference": sorted(cpp), "java": sorted(java),
                 "missing_in_java": sorted(cpp - java), "extra_in_java": sorted(java - cpp)}
        result["options"][tool] = entry
        result["gaps"] += len(entry["missing_in_java"]) + len(entry["extra_in_java"])

    declared, defined = {}, {}
    for unit, (headers, sources, _) in UNITS.items():
        for header in headers:
            unit_declared, unit_defined = cpp_methods(SOURCE / header, [SOURCE / s for s in sources],
                                                      Path(header).stem)
            declared[unit] = declared.get(unit, set()) | unit_declared
            defined[unit] = defined.get(unit, set()) | unit_defined
    every = set().union(*defined.values())
    static = set().union(*(defined[unit] for unit in STATIC_UNITS))
    live = reachable(every, static, call_graph())

    for unit, (headers, sources, java_files) in UNITS.items():
        ported = java_methods(java_files)
        missing = sorted(defined[unit] - ported)
        unreachable_here = sorted(name for name in missing if name not in live)
        renamed, unexplained = {}, []
        for name in missing:
            if name in unreachable_here:
                continue
            if name in PORTED_AS:
                file, symbol = PORTED_AS[name]
                if symbol in java_methods([file]):
                    renamed[name] = "%s: %s()" % (file, symbol)
                    continue
            unexplained.append(name)
        entry = {
            "reference_methods": len(defined[unit]),
            "ported_by_name": len(defined[unit] & ported),
            "ported_under_another_name": renamed,
            "defined_but_unreachable_in_reference": unreachable_here,
            "missing_in_java": unexplained,
            "declared_but_never_defined_in_reference": sorted(declared[unit] - defined[unit]),
            "java_only": sorted(ported - defined[unit]),
        }
        result["methods"][unit] = entry
        result["ported_as"].update(renamed)
        result["gaps"] += len(unexplained)

    unreachable = {name for entry in result["methods"].values()
                   for name in entry["defined_but_unreachable_in_reference"]}
    result["reachable"] = len(every - unreachable)
    result["ported"] = result["reachable"] - result["gaps"]

    if not options.no_write:
        OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    print("OPTIONS", {tool: "%d/%d" % (len(entry["reference"]) - len(entry["missing_in_java"]),
                                       len(entry["reference"])) for tool, entry in result["options"].items()})
    for unit, entry in result["methods"].items():
        print("METHODS %-24s %3d by name, %d renamed/inlined, %d unreachable in reference, missing=%s" % (
            unit, entry["ported_by_name"], len(entry["ported_under_another_name"]),
            len(entry["defined_but_unreachable_in_reference"]), entry["missing_in_java"] or "[]"))
    if options.print:
        print(json.dumps(result, indent=2))
    print("GAPS %d" % result["gaps"])
    print("REPORT " + OUTPUT.relative_to(ROOT).as_posix())
    return 1 if result["gaps"] else 0


if __name__ == "__main__":
    sys.exit(main())
