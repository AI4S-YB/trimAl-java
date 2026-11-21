#!/bin/bash
#
# SPDX-License-Identifier: GPL-3.0-or-later
# This file is part of the trimAl Java port verification tooling
# (https://github.com/inab/trimal, GPL-3.0-or-later).
#
# Parallel ranged downloader for servers that throttle per connection (~10KB/s each).
# Usage: pdl.sh <url> <size> <outfile> <parts> <workdir>
set -u
URL="$1"; SIZE="$2"; OUTFILE="$3"; N="$4"; WORK="$5"
mkdir -p "$WORK"
CHUNK=$(( (SIZE + N - 1) / N ))

dl_chunk() {
  i=$1; S=$((i * CHUNK)); E=$((S + CHUNK - 1))
  if [ "$E" -ge "$SIZE" ]; then E=$((SIZE - 1)); fi
  OUT="$WORK/part_$(printf %03d $i)"
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

TOTAL=0
for i in $(seq 0 $((N - 1))); do
  SZ=$(stat -c %s "$WORK/part_$(printf %03d $i)" 2>/dev/null || echo 0)
  TOTAL=$((TOTAL + SZ))
done
echo "downloaded=$TOTAL / $SIZE"
if [ "$TOTAL" -eq "$SIZE" ]; then
  cat "$WORK"/part_* > "$OUTFILE"
  echo "ASSEMBLED $OUTFILE"
fi
