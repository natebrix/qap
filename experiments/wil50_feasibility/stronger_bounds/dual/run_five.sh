#!/bin/bash
# research_bound_quality.md candidate 2: the regional 5-point count inequalities
# (latest/dnncut.h, dnncuts -five K) against the triangles alone, same protocol
# both arms: from the give-up point (400 iterations), 12 rounds of 300 with
# 1,000 triangles a round, Polyak toward best + 50 uncapped, then a quiet round
# of 1,500 iterations from the best round's multipliers. The five arm adds 500
# of the 5-point cuts a round. Nodes (all32 line numbers; d<k> = depth4/
# nodes.txt row k): 14 (84 short at convergence, 8.7 short after theta 0.4 +
# quiet), 10 (158 short), d2 (109 short), 13 (182 short). 4 at a time, one
# core each; resumable (FINAL line = done). Logs: {tri,five}_<id>.log.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
C="$R/experiments/dnn_cuts/dnncuts -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 -i 300 -r 12 -delta 50 -nocap -freeze 0:1500"
spec() { case $1 in d*) awk -v k=${1#d} '$1==k{print $3}' "$D/../depth4/nodes.txt";; *) sed -n ${1}p "$D/../../tuning/giveup/all32.txt";; esac; }
one() {
  name=$1; flags=${2//_/ }; [ "$flags" = " " ] && flags=""; k=$3; spec=$4
  f=$D/${name}_$k.log
  grep -q '^FINAL' "$f" 2>/dev/null && exit 0
  $C -x "$spec" $flags "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
export -f one; export D R C
{
  for k in 14 10 d2 13; do
    echo "tri - $k $(spec $k)"
    echo "five -five_500 $k $(spec $k)"
  done
} | sed 's/^tri - /tri _ /' | xargs -P 4 -L 1 bash -c 'one "$@"' _
