#!/usr/bin/env bash
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#

# Build with unit tests, then run the verification gates through verify/run_full.py.
# Everything after the build is the driver's: `--only quick,numeric`, `--skip pffull`,
# `--skip-missing` and its other options are passed on. Without arguments every gate
# runs, including the slow ones and those that need downloaded data or a Windows JVM.
#
# Launched from Git Bash, the JAR is built on the Windows side and the gates run
# under WSL, where the reference executables work. No Maven installation is needed:
# the build goes through the Maven wrapper.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

if command -v cygpath >/dev/null 2>&1 && [[ "${TRIMAL_SKIP_BUILD:-0}" != 1 ]]; then
  java/mvnw -B -f java/pom.xml package
  WSL_ROOT="$(wsl -e wslpath -a "$ROOT" </dev/null | tr -d '\r')"
  TRIMAL_SKIP_BUILD=1 exec wsl -e bash "$WSL_ROOT/verify/run_all.sh" "$@" </dev/null
fi

if [[ "${TRIMAL_SKIP_BUILD:-0}" != 1 ]]; then
  java/mvnw -B -f java/pom.xml package
fi
# The repository checks: every test resource is tracked, scripts carry their
# executable bit, no tool names one machine, and the workflow files are valid. An
# invalid workflow fails every CI run in zero seconds, with no job and no log, so
# this is the only place that can catch it before a push.
python3 verify/repo_check.py
python3 verify/run_full.py "$@"
