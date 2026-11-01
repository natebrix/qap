#!/bin/bash
# Follow-up: both fixes together (uncapped target + quiet round) on the
# depth-4 children 1, 4, 11 (62-66 short at convergence; with the quiet
# round alone child 1 ended 11.5 short) and on node 14 (84 short; uncapped
# alone ended 5 short). Same protocol otherwise (give-up point, 300-iteration
# rounds, delta 50, 12 rounds). Logs: both_<id>.log; resumable.
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
for k in d1 d4 d11 14; do echo "both -I_400_-i_300_-r_12_-delta_50_-nocap_-freeze_0:600 $k $(spec $k)"; done | xargs -P 4 -L 1 bash -c 'one "$@"' _
