#!/bin/bash
# Cost experiment: the four near misses the regional cuts fathom (step 2:
# nodes 6, 9, 15, 16, fathomed in rounds 5-7 of 700 iterations, 5,000-6,400
# iterations and 2,850-3,740 CPU-s each) under cheaper protocols. Question:
# how few iterations does the cut phase need, from a converged uncut bound
# and from the give-up point the search would hand over (about 400
# iterations for these nodes under H = 300).
#   i300     -I 1500 -i 300 -r 10 -delta 100   fewer iterations per round
#   i300d30  -I 1500 -i 300 -r 10 -delta 30    ... and a damped step
#   g400     -I  400 -i 300 -r 12 -delta 50    from the give-up point
#   g400s    -I  400 -i 150 -r 20 -delta 50    ... with short rounds
# 4 at a time, one core each; dnncuts stops at the first certified bound
# above the cutoff. Logs: <proto>_<id>.log; resumable (FINAL line = done).
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
PROTOS=${PROTOS:-"i300:-I_1500_-i_300_-r_10_-delta_100 i300d30:-I_1500_-i_300_-r_10_-delta_30 g400:-I_400_-i_300_-r_12_-delta_50 g400s:-I_400_-i_150_-r_20_-delta_50"}
NODES=${NODES:-"6 9 15 16"}
one() {
  name=$1; flags=${2//_/ }; k=$3; spec=$4
  f=$D/${name}_$k.log
  grep -q '^FINAL' "$f" 2>/dev/null && exit 0
  "$R/experiments/dnn_cuts/dnncuts" -grid 5:10 -x "$spec" -u 48816 -T 1 -k 1000 $flags -v \
    "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
export -f one; export D R
for p in $PROTOS; do
  for k in $NODES; do echo "${p%%:*} ${p#*:} $k $(sed -n ${k}p "$D/../../tuning/giveup/all32.txt")"; done
done | xargs -P 4 -L 1 bash -c 'one "$@"' _
