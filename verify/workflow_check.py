#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Check the workflow files before they are pushed.

GitHub validates a workflow file before it starts anything. An invalid file fails
every run in zero seconds with no job and no log, on every branch at once, and
nothing in CI can report it: a typo in an expression once did exactly that. This
script parses each workflow as YAML and checks what GitHub's first pass checks:
the structure of `on`, `jobs`, `steps`, matrix references, and that every
`${{ ... }}` is closed before the next one opens and names known contexts.

    python3 verify/workflow_check.py
"""
from pathlib import Path
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
WORKFLOWS = ROOT / ".github" / "workflows"
EXPRESSION = re.compile(r"\$\{\{(.*?)\}\}", re.S)
CONTEXTS = {"github", "env", "job", "jobs", "steps", "runner", "secrets", "inputs", "matrix",
            "needs", "vars"}
STEP_KEYS = {"name", "id", "if", "uses", "run", "shell", "working-directory", "env", "with",
             "continue-on-error", "timeout-minutes", "retry"}
# expressions GitHub evaluates while it still is reading the file: they may only use
# the contexts that exist at that point
STATIC_CONTEXTS = {"github", "inputs", "vars"}


def walk(node, path, problems):
    """Every string in the workflow: closed expressions, known contexts."""
    if isinstance(node, dict):
        for key, value in node.items():
            walk(value, path + [str(key)], problems)
    elif isinstance(node, list):
        for index, value in enumerate(node):
            walk(value, path + ["[%d]" % index], problems)
    elif isinstance(node, str):
        opened = node.count("${{")
        if opened != len(EXPRESSION.findall(node)):
            problems.append("%s: an expression is not closed: %r" % (".".join(path), node))
        for expression in EXPRESSION.findall(node):
            name = expression.strip().split(".")[0].split("[")[0].strip()
            if name and name not in CONTEXTS and not name.isdigit() \
                    and name not in ("true", "false", "null") and not name.startswith(("'", '"')):
                problems.append("%s: unknown context %r in %r" % (".".join(path), name, expression.strip()))


def check(path):
    problems = []
    try:
        import yaml
    except ImportError:
        yaml = None
    text = path.read_text(encoding="utf-8")
    if yaml is None:
        for number, line in enumerate(text.splitlines(), 1):
            if line.count("${{") != len(EXPRESSION.findall(line)):
                problems.append("%s:%d: an expression is not closed: %r" % (path.name, number, line))
        return problems
    try:
        document = yaml.safe_load(text)
    except yaml.YAMLError as error:
        return ["%s: not valid YAML: %s" % (path.name, " ".join(str(error).split()))]
    if not isinstance(document, dict):
        return ["%s: not a mapping" % path.name]
    # YAML 1.1 reads the key `on` as the boolean true
    on = document.get(True, document.get("on"))
    if on is None:
        problems.append("%s: no `on` section" % path.name)
    elif not isinstance(on, (dict, list, str)):
        problems.append("%s: `on` must be a mapping, a list or a string" % path.name)
    elif isinstance(on, dict):
        for event, value in on.items():
            if value is not None and not isinstance(value, (dict, list)):
                problems.append("%s: on.%s must be empty, a mapping or a list" % (path.name, event))
    jobs = document.get("jobs")
    if not isinstance(jobs, dict) or not jobs:
        problems.append("%s: no jobs" % path.name)
        return problems
    for name, job in jobs.items():
        where = "%s: jobs.%s" % (path.name, name)
        if not isinstance(job, dict):
            problems.append(where + ": not a mapping")
            continue
        if "uses" in job:
            continue                      # a reusable workflow has no steps of its own
        if "runs-on" not in job:
            problems.append(where + ": no runs-on")
        steps = job.get("steps")
        if not isinstance(steps, list) or not steps:
            problems.append(where + ": no steps")
            continue
        matrix = (job.get("strategy") or {}).get("matrix") or {}
        for index, step in enumerate(steps):
            place = "%s: steps[%d]" % (where, index)
            if not isinstance(step, dict):
                problems.append(place + ": not a mapping")
                continue
            for key in step:
                if key not in STEP_KEYS:
                    problems.append("%s: unknown key %r" % (place, key))
            if "uses" in step and "run" in step:
                problems.append(place + ": a step has either uses or run")
            if "uses" not in step and "run" not in step:
                problems.append(place + ": a step has neither uses nor run")
            if "uses" in step and not re.fullmatch(r"[\w.-]+/[\w.-]+@[\w./-]+", str(step["uses"])):
                problems.append("%s: malformed uses %r" % (place, step["uses"]))
            for key in ("if", "name", "run", "working-directory"):
                if key in step:
                    walk(step[key], [place, key], problems)
            walk(step.get("with", {}), [place, "with"], problems)
            walk(step.get("env", {}), [place, "env"], problems)
            for expression in EXPRESSION.findall(json.dumps(step, default=str)):
                if "matrix." in expression and not matrix:
                    problems.append("%s: uses matrix but the job has none" % place)
                elif "matrix." in expression:
                    for axis in set(re.findall(r"matrix\.([\w-]+)", expression)):
                        if axis not in matrix and axis != "github":
                            problems.append("%s: matrix has no axis %r" % (place, axis))
    # expressions in the file-level sections are evaluated while the file is read
    for section in ("name", "concurrency", "run-name"):
        if section in document:
            walk(document[section], ["%s:%s" % (path.name, section)], problems)
    walk(document.get("env", {}), ["%s:env" % path.name], problems)
    for job in jobs.values():
        if not isinstance(job, dict):
            continue
        matrix = (job.get("strategy") or {}).get("matrix") or {}
        for key in ("name", "runs-on", "if", "timeout-minutes", "continue-on-error"):
            if key in job:
                walk(job[key], ["%s:jobs.*.%s" % (path.name, key)], problems)
                for expression in EXPRESSION.findall(json.dumps(job[key], default=str)):
                    for axis in set(re.findall(r"matrix\.([\w-]+)", expression)):
                        if axis not in matrix and axis != "github":
                            problems.append("%s: jobs.*.%s: matrix has no axis %r" % (path.name, key, axis))
    return problems


def main():
    problems = []
    for path in sorted(WORKFLOWS.glob("*.yml")):
        problems += check(path)
    for problem in problems:
        print("PROBLEM", problem)
    print("WORKFLOWS %d file(s), %d problem(s)" % (len(list(WORKFLOWS.glob("*.yml"))), len(problems)))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
