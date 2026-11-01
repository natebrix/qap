#!/bin/bash
# First experiments of research_bound_quality.md section 3.2: the cheap fixes
# to the dual method, on the nodes of all32.txt (ids = line numbers) and the
# depth-4 children of depth4/nodes.txt (d<id>). 4 at a time, one core each.
#   nocap    near misses 6, 9, 15, 16 and node 14 from the give-up point
#            (step-3 protocol) with the Polyak target best + delta instead of
#            min(cutoff, best + delta): fathom round and CPU against step 3
#   quiet150 the 150-iteration protocol that failed in step 3, plus a quiet
#            round of 600 iterations (multipliers frozen, no new cuts)
#   quiet4   depth-4 children 1, 2, 4, 11 (62-109 short) from the give-up
#            point, 300-iteration rounds, plus a quiet round of 600: do they
#            fathom once the certificate catches up with the Lagrangian value?
#   k5000    clear misses 2 and 7, step-2 protocol with 5,000 cuts per round
#   t04      overshooters 3 and 8, step-2 protocol with theta 0.4
# Logs: <proto>_<id>.log; resumable (FINAL line = done).
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
{
  for k in 6 9 15 16 14; do echo "nocap -I_400_-i_300_-r_12_-delta_50_-nocap $k $(spec $k)"; done
  for k in 6 9 15 16; do echo "quiet150 -I_400_-i_150_-r_20_-delta_50_-freeze_0:600 $k $(spec $k)"; done
  for k in d1 d2 d4 d11; do echo "quiet4 -I_400_-i_300_-r_12_-delta_50_-freeze_0:600 $k $(spec $k)"; done
  for k in 2 7; do echo "k5000 -I_1500_-i_700_-r_8_-delta_100_-k_5000 $k $(spec $k)"; done
  for k in 3 8; do echo "t04 -I_1500_-i_700_-r_8_-delta_100_-t_0.4 $k $(spec $k)"; done
} | xargs -P 4 -L 1 bash -c 'one "$@"' _
