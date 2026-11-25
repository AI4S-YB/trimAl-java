#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Repository checks that a Windows workstation cannot make by itself, and CI runs
first.

    python3 verify/repo_check.py

Three things went wrong exactly once each before this existed:

- an unanchored `upstream/` in .gitignore silently kept the golden inputs out of
  Git: they were on the developer's disk, so the tests passed there and failed on
  every clean checkout;
- on Windows `core.filemode` is off, so no script ever received its executable
  bit, and `./mvnw` was refused on macOS;
- an invalid workflow file fails every CI run in zero seconds with no log.

So this checks, in the index (not in the working tree, which is what the developer
sees):

1. every file under java/src/test/resources/ that exists on disk is tracked;
2. the launchers, shell scripts and the Maven wrapper carry mode 100755;
3. no text file is stored with CRLF line endings;
4. no tracked tool file names a machine of its own (/mnt/<drive>, a home directory,
   a checked-out path);
5. the workflow files are valid (verify/workflow_check.py).
"""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "verify"))

import workflow_check  # noqa: E402

RESOURCES = ROOT / "java" / "src" / "test" / "resources"
EXECUTABLES = ["java/mvnw", "java/trimal", "java/readal", "java/statal",
               "verify/run_all.sh", "verify/compare_trimal.sh", "verify/pdl.sh"]
# tracked tool files that must not name one machine
TOOL_PATHS = ["verify/*.py", "verify/*/*.py", "verify/*/*/*.py", "verify/*.sh", "java/*.cmd",
              "java/src/**/*.java", ".github/**/*.yml", "java/pom.xml"]


def git(*arguments):
    return subprocess.run(["git"] + list(arguments), cwd=ROOT, capture_output=True, text=True,
                          encoding="utf-8", errors="replace")


def resources_tracked(problems):
    if not RESOURCES.is_dir():
        return
    tracked = set(git("ls-files", "java/src/test/resources").stdout.splitlines())
    for path in sorted(p for p in RESOURCES.rglob("*") if p.is_file()):
        relative = path.relative_to(ROOT).as_posix()
        if relative not in tracked:
            ignored = git("check-ignore", "-v", relative).stdout.strip()
            problems.append("test resource is not tracked: %s%s" % (
                relative, " (%s)" % ignored.split("\t")[0] if ignored else ""))


def executable_bits(problems):
    listing = git("ls-files", "-s", "--", *EXECUTABLES)
    modes = dict((line.split("\t", 1)[1].strip(), line.split(" ", 1)[0])
                 for line in listing.stdout.splitlines() if "\t" in line)
    for relative in EXECUTABLES:
        if relative not in modes:
            problems.append("not tracked: " + relative)
        elif modes[relative] != "100755":
            problems.append("mode %s, expected 100755: %s" % (modes[relative], relative))


def stored_line_endings(problems):
    """Text files are stored with LF, whatever a checkout makes of them (.gitattributes).
    A commit that stores CRLF shows up as modified in every clean checkout, which is
    how CI's "runners left every tracked file untouched" step once failed on a
    dependency update. Files marked -text (byte-exact test fixtures) are exempt."""
    for line in git("ls-files", "--eol").stdout.splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[0] in ("i/crlf", "i/mixed") and "-text" not in fields[2]:
            problems.append("stored with CRLF, expected LF: " + fields[-1])


def machine_markers(problems):
    # this file spells out the patterns it searches for
    found = git("grep", "-n", "-I", "-e", "ClaudeCodeDev", "-e", "/mnt/[a-z]/", "-e", "/home/",
                "-e", "[A-Z]:\\\\", "-e", "Program Files", "--", *TOOL_PATHS,
                ":(exclude)verify/repo_check.py")
    for line in found.stdout.splitlines():
        problems.append("names a machine: " + line[:160])


def main():
    problems = []
    resources_tracked(problems)
    executable_bits(problems)
    stored_line_endings(problems)
    machine_markers(problems)
    for path in sorted((ROOT / ".github" / "workflows").glob("*.yml")):
        problems += workflow_check.check(path)
    for problem in problems:
        print("PROBLEM", problem)
    print("REPO %d problem(s)" % len(problems))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
