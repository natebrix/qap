#!/bin/bash
# Give-up threshold A/B on the cloud (Xeon): the 32 nodes of all32.txt with
# -H 300 and -H 200, solver mode with PRSM and degree 4, -i 2500 -u 48816,
# 4 at a time. Resumable: a node whose log already has a RESULT line is skipped.
# Logs: cloud/H<H>_<id>.log, id = line number in all32.txt
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
one() {
  H=$1; k=$2; spec=$3
  f=$D/cloud/H${H}_$k.log
  grep -q '^RESULT' "$f" 2>/dev/null && exit 0
  "$R/latest/dnntest" -T 1 -d 4 -g 1.0 -G 1.0 -H $H -i 2500 -u 48816 -x "$spec" "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
export -f one; export D R
for H in 300 200; do awk -v H=$H '{print H, NR, $1}' "$D/all32.txt"; done |
  xargs -P 4 -L 1 bash -c 'one "$@"' _
