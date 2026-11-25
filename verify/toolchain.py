#!/usr/bin/env python3
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

"""Locate the build and run tools the verification runners call.

The runners used to call /usr/bin/java and /usr/bin/javac directly, which bypasses
the JDK a CI setup step selects (and whatever a developer put first on PATH).
java_tool() prefers $JAVA_HOME/bin and falls back to PATH; tool() resolves every
other program on PATH and fails with a clear message instead of assuming /usr/bin.
"""
import functools
import os
from pathlib import Path
import shutil


@functools.lru_cache(maxsize=None)
def tool(name):
    found = shutil.which(name)
    if found is None:
        raise SystemExit("required tool not found on PATH: " + name)
    return found


@functools.lru_cache(maxsize=None)
def java_tool(name):
    home = os.environ.get("JAVA_HOME")
    if home:
        candidate = Path(home) / "bin" / name
        if candidate.is_file():
            return str(candidate)
    return tool(name)
