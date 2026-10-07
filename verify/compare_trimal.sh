#!/usr/bin/env bash
# Run the strict native/Java differential test on Linux or through WSL from Git Bash.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="${1:-verify/quick_manifest.txt}"
if [[ $# -gt 0 ]]; then shift; fi
if command -v cygpath >/dev/null 2>&1; then
  WINDOWS_ROOT="$(cygpath -m "$ROOT")"
  WSL_ROOT="$(wsl -e wslpath -a "$WINDOWS_ROOT" </dev/null | tr -d '\r')"
  MSYS_NO_PATHCONV=1 wsl -e bash -c 'cd "$1"; shift; exec python3 verify/differential.py "$@"' _ "$WSL_ROOT" "$MANIFEST" "$@" </dev/null
else
  cd "$ROOT"
  exec python3 verify/differential.py "$MANIFEST" "$@" </dev/null
fi
