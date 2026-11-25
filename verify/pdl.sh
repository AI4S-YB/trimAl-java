#!/bin/bash
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#
# Parallel ranged downloader, for servers that throttle each connection. The range is
# split into <parts> chunks that are fetched concurrently and resumed where they
# stopped; the chunks are then assembled into <outfile>, whose size is checked.
# Usage: pdl.sh <url> <size-in-bytes> <outfile> <parts> <workdir>
set -u
URL="$1"; SIZE="$2"; OUTFILE="$3"; N="$4"; WORK="$5"
mkdir -p "$WORK"
CHUNK=$(( (SIZE + N - 1) / N ))

dl_chunk() {
  i=$1; S=$((i * CHUNK)); E=$((S + CHUNK - 1))
  if [ "$E" -ge "$SIZE" ]; then E=$((SIZE - 1)); fi
  OUT="$WORK/part_$(printf %03d "$i")"
  EXPECT=$((E - S + 1))
  for try in 1 2 3 4 5 6 7 8 9 10; do
    SZ=$(stat -c %s "$OUT" 2>/dev/null || echo 0)
    if [ "$SZ" -ge "$EXPECT" ]; then return; fi
    curl -sS --max-time 3600 -r $((S + SZ))-$E -o - "$URL" >> "$OUT" 2>/dev/null
    sleep 1
  done
}

for i in $(seq 0 $((N - 1))); do
  dl_chunk "$i" &
done
wait

total=0
for i in $(seq 0 $((N - 1))); do
  part="$WORK/part_$(printf %03d "$i")"
  total=$((total + $(stat -c %s "$part" 2>/dev/null || echo 0)))
done
if [ "$total" -ne "$SIZE" ]; then
  echo "incomplete download: $total of $SIZE bytes under $WORK" >&2
  exit 1
fi
cat "$WORK"/part_* > "$OUTFILE"
echo "ASSEMBLED $OUTFILE ($SIZE bytes)"
