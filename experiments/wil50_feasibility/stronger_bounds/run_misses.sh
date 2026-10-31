#!/bin/bash
# Step 2 of the regional-cut evaluation: the 14 misses of
# ../tuning/giveup/all32.txt (the nodes degree-4 PRSM does not fathom under
# the give-up rule; ids = line numbers), each run through dnncuts with the
# grid regions of wil50, 4 at a time. The uncut reference is round 0 of the
# same log (first_lb, the bound after -I iterations) and, for the equal-budget
# question, the control runs of node 1 (node1/control.log).
# Configurations (name:flags, '_' for spaces), set CONFIGS to override:
#   polyak   8 rounds of 700 iterations, Polyak step toward best + 100, new
#            cuts at multiplier 0 (the node-1 protocol)
# Resumable: a node whose log has a FINAL line is skipped.
# Logs: misses/<name>_<id>.log
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../.." && pwd); mkdir -p "$D/misses"
CONFIGS=${CONFIGS:-"polyak:-I_1500_-i_700_-r_8_-k_1000_-delta_100"}
MISSES=${MISSES:-"1 2 3 4 5 6 7 8 9 10 13 14 15 16"}
one() {
  name=$1; flags=${2//_/ }; k=$3; spec=$4
  f=$D/misses/${name}_$k.log
  grep -q '^FINAL' "$f" 2>/dev/null && exit 0
  "$R/experiments/dnn_cuts/dnncuts" -grid 5:10 -x "$spec" -u 48816 -T 1 $flags -v \
    "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
export -f one; export D R
for c in $CONFIGS; do
  for k in $MISSES; do echo "${c%%:*} ${c#*:} $k $(sed -n ${k}p "$D/../tuning/giveup/all32.txt")"; done
done | xargs -P 4 -L 1 bash -c 'one "$@"' _
