#!/bin/bash
# The local-7 inequalities as a pool cut type (latest/dnncut.C type 6, dnncuts
# -seven K[:S]): the step-7 protocol on nodes 13 and d2 with the local-7 cuts
# on top of the regional triangles, and alone. References: ../dual/tri_13.log
# (48,755.42 after the quiet round, 58.6 short), ../dual/tri_d2.log (48,789.87,
# 24.2 short), both 12 rounds of 1,000 triangles, uncapped, quiet 1,500 from
# the best round's multipliers.
#   triseven   1,000 triangles + 300 local-7 cuts a round (1,000 random
#              7-subsets per region and round, plus a third as many from the
#              15 most fractional elements)
#   triseven3k the same with 3,000 subsets per region (more separation)
#   sevenonly  300 local-7 cuts a round, no triangles
# One core each, 2 at a time. Logs <proto>_<id>.log; FINAL line = done.
#   bash run_pool.sh [13 d2]
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
C="$R/experiments/dnn_cuts/dnncuts -grid 5:10 -u 48816 -T 1 -v -I 400 -i 300 -r 12 -delta 50 -nocap -freeze 0:1500"
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
    echo "triseven -k_1000_-seven_300:1000 $k"
    echo "triseven3k -k_1000_-seven_300:3000 $k"
    echo "sevenonly -k_0_-seven_300:1000 $k"
  done
} | xargs -P 2 -L 1 bash -c 'one "$@"' _
for f in $D/*_1[03].log $D/*_d2.log; do [ -f "$f" ] && echo "$(basename $f): $(grep '^FINAL' $f) | quiet: $(grep '^quiet' $f | grep -o 'lb [0-9.]*') | rounds cpu: $(grep '^round' $f | grep -o 'cpu [0-9.]*' | tail -1)"; done
