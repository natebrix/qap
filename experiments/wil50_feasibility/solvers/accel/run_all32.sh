#!/bin/bash
# A/B of the iteration accelerations on the 32 wil50 nodes of
# tuning/giveup/all32.txt: solver mode with PRSM (-g 1.0 -G 1.0), -H 300,
# -i 2500 -u 48816, 4 at a time. Configurations (name:flags):
#   d4     degree 4, Weyl lam_hi (the baseline)
#   d3     degree 3
#   d3L2   degree 3, lam_hi from 2 power steps (-L 2)
#   d4L2   degree 4, -L 2
#   d4W5   degree 4, Anderson acceleration, memory 5 (-W 5)
#   d4A25  degree 4, Halpern anchoring, restarts at most 25 apart (-A 25)
# Resumable: a node whose log has a RESULT line is skipped.
# Logs: runs/<name>_<id>.log, id = line number in all32.txt
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
CONFIGS=${CONFIGS:-"d4:-d_4 d3:-d_3 d3L2:-d_3_-L_2 d4L2:-d_4_-L_2 d4W5:-d_4_-W_5 d4A25:-d_4_-A_25"}
one() {
  name=$1; flags=${2//_/ }; k=$3; spec=$4
  f=$D/runs/${name}_$k.log
  grep -q '^RESULT' "$f" 2>/dev/null && exit 0
  "$R/latest/dnntest" -T 1 $flags -g 1.0 -G 1.0 -H 300 -i 2500 -u 48816 -v -x "$spec" \
    "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
export -f one; export D R
for c in $CONFIGS; do
  awk -v n=${c%%:*} -v fl=${c#*:} '{print n, fl, NR, $1}' "$D/../../tuning/giveup/all32.txt"
done | xargs -P 4 -L 1 bash -c 'one "$@"' _
