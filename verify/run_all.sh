#!/usr/bin/env bash
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#
# Build with tests, generate the quick matrix, and compare all output bytes.
# Pass --scale to additionally run the slow 3,583-sequence conservation jobs.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if command -v cygpath >/dev/null 2>&1; then
  cd "$ROOT"
  mvn -B -f java/pom.xml package
  WINDOWS_ROOT="$(cygpath -m "$ROOT")"
  WSL_ROOT="$(wsl -e wslpath -a "$WINDOWS_ROOT" </dev/null | tr -d '\r')"
  MSYS_NO_PATHCONV=1 exec wsl -e env TRIMAL_SKIP_BUILD=1 bash "$WSL_ROOT/verify/run_all.sh" "$@" </dev/null
fi
cd "$ROOT"
if [[ "${TRIMAL_SKIP_BUILD:-0}" != 1 ]]; then
  mvn -B -f java/pom.xml package
fi
python3 verify/generate_suite.py
python3 verify/differential.py verify/quick_manifest.txt --timeout 90 --jobs 2
python3 verify/differential.py verify/extended_manifest.txt --timeout 90 --jobs 2
python3 verify/differential.py verify/legacy_manifest.txt --timeout 90 --jobs 2
python3 verify/companions/verify_companions.py --java-jar "$ROOT/java/target/trimal.jar"
if [[ "${1:-}" == --scale ]]; then
  python3 verify/differential.py verify/scale_manifest.txt --timeout 3600 --jobs 1
fi
