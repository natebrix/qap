#!/bin/bash
# The dual-method test (research_bound_quality.md section 4, protocol B): how
# far below the cut pool's dual optimum do the 12 Polyak rounds of step 7 stop?
# Nodes 13 (21:28,8:41,28:19; tri alone +158.6, 58.6 short after the quiet
# round, dual/tri_13.log) and d2 (39:27,30:48,25:49,26:5; +133.0, 24.2 short,
# dual/tri_d2.log). Same start as step 7 (400 uncut iterations, 1,000 regional
# triangles a round, uncapped Polyak, quiet round of 1,500 from the best
# round's multipliers), three protocols:
#   slow30        30 rounds of 300 with a small step (theta 0.3, target best+30)
#   slow30_pool12 the same, no new cuts after round 12 (-addto 12): the gain
#                 from better multipliers on the step-7 pool alone
#   fast60        60 rounds of 100 with 500 cuts a round, same total budget:
#                 frequent small updates, a proxy for multipliers moved inside
#                 the ADMM
# One core each; run 2 at a time. Logs <proto>_<id>.log; FINAL line = done.
#   bash run_dualtest.sh            (all six)
#   bash run_dualtest.sh 13         (one node)
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
C="$R/experiments/dnn_cuts/dnncuts -grid 5:10 -u 48816 -T 1 -v -I 400 -nocap -freeze 0:1500"
spec() { case $1 in d2) echo 39:27,30:48,25:49,26:5;; 13) echo 21:28,8:41,28:19;; 14) echo 47:23,10:44,40:34;; 10) echo 39:27,30:48,25:49;; esac; }
one() {
  name=$1; flags=${2//_/ }; k=$3
  f=$D/${name}_$k.log
  grep -q '^FINAL' "$f" 2>/dev/null && exit 0
  $C -x "$(spec $k)" $flags "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
export -f one spec; export D R C
nodes=${@:-13 d2}
{
  for k in $nodes; do
    echo "slow30 -i_300_-r_30_-k_1000_-delta_30_-t_0.3 $k"
    echo "slow30_pool12 -i_300_-r_30_-k_1000_-delta_30_-t_0.3_-addto_12 $k"
    echo "fast60 -i_100_-r_60_-k_500_-delta_30_-t_0.3 $k"
  done
} | xargs -P 2 -L 1 bash -c 'one "$@"' _
for f in $D/*_*.log; do echo "$(basename $f): $(grep '^FINAL' $f) | quiet: $(grep '^quiet' $f | grep -o 'lb [0-9.]*')"; done
