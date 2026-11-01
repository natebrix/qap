#!/bin/bash
# Depth-4 children of the clear misses (nodes.txt from make_sample.py), two
# runs each, 4 at a time:
#   uncut_<id>.log  dnntest, PRSM + degree 4, no give-up, 1,500 iterations
#                   (-v: the certified bound every 10 iterations, so the
#                   give-up-point value at 400 is in the trace)
#   cuts_<id>.log   dnncuts from the give-up point: 400 uncut iterations, then
#                   300-iteration rounds with 1,000 regional triangles each,
#                   Polyak step toward best + 50, up to 12 rounds, stopping at
#                   the first certified bound above the cutoff
# Resumable (RESULT / FINAL line = done).
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
W=$R/qaplib/wil50.dat
one() {
  kind=$1; k=$2; spec=$3
  f=$D/${kind}_$k.log
  grep -q '^RESULT\|^FINAL' "$f" 2>/dev/null && exit 0
  if [ "$kind" = uncut ]; then
    "$R/latest/dnntest" -T 1 -d 4 -g 1.0 -G 1.0 -H 0 -i 1500 -u 48816 -v -x "$spec" "$W" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
  else
    "$R/experiments/dnn_cuts/dnncuts" -grid 5:10 -x "$spec" -u 48816 -T 1 -k 1000 -I 400 -i 300 -r 12 -delta 50 -v "$W" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
  fi
}
export -f one; export D R W
{ awk '{print "uncut", $1, $3}' "$D/nodes.txt"; awk '{print "cuts", $1, $3}' "$D/nodes.txt"; } | xargs -P 4 -L 1 bash -c 'one "$@"' _
