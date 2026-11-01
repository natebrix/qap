#!/bin/bash
# theta 0.4 + quiet round on the two nodes the quiet round left 8-12 short
# (depth-4 child 1, 62 short at convergence; node 14, 84 short): the
# overshoot rounds cost them 2-3 rounds each. Logs: tq_<id>.log; resumable.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
C="$R/experiments/dnn_cuts/dnncuts -grid 5:10 -u 48816 -T 1 -k 1000 -v"
spec() { case $1 in d*) awk -v k=${1#d} '$1==k{print $3}' "$D/../depth4/nodes.txt";; *) sed -n ${1}p "$D/../../tuning/giveup/all32.txt";; esac; }
one() {
  name=$1; flags=${2//_/ }; k=$3; spec=$4
  f=$D/${name}_$k.log
  grep -q '^FINAL' "$f" 2>/dev/null && exit 0
  $C -x "$spec" $flags "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
export -f one; export D R C
for k in d1 14; do echo "tq -I_400_-i_300_-r_12_-delta_50_-t_0.4_-nocap_-freeze_0:600 $k $(spec $k)"; done | xargs -P 4 -L 1 bash -c 'one "$@"' _
