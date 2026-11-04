#!/bin/bash
# The box (count) inequalities as an in-pool cut type (dnncut.h, cut type 8;
# new_cuts/literature/REPORT.md, recommendation 2). Predeclared protocol,
# step 7's (stronger_bounds/dual/run_five.sh) with the search's per-node
# settings (bundle D): from the give-up point (400 iterations), 12 Polyak
# rounds of 300 iterations with 1,000 triangles a round, target best + 50
# uncapped, quiet round of 1,500 from the best round's multipliers. The box
# arm adds 200 box cuts a round (|Q| <= 24, 40 starts per region and b).
# Nodes 10, 13, d2 (the gate) and 14 (fathomed by the triangles in step 7: a
# check). Gate: KILL if the median over 10, 13, d2 of (box - tri) retained
# certificates is below +15. Each node's two arms run side by side, so each
# pair shares a host; the host's CPU model is logged. Resumable (FINAL = done).
#   bash run_box_pool.sh            # both waves, 4 processes at a time
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
export DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
C="$R/experiments/dnn_cuts/dnncuts -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 -i 300 -r 12 -delta 50 -nocap -freeze 0:1500"
spec() { case $1 in d*) awk -v k=${1#d} '$1==k{print $3}' "$R/experiments/wil50_feasibility/stronger_bounds/depth4/nodes.txt";;
                    *) sed -n ${1}p "$R/experiments/wil50_feasibility/tuning/giveup/all32.txt";; esac; }
one() {  # arm node core
  arm=$1; k=$2; core=$3; f=$D/runs/${arm}_$k.log
  grep -q '^FINAL' "$f" 2>/dev/null && return
  extra=""; [ "$arm" = box ] && extra="-box 200"
  { echo "host $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2) amx $(grep -c amx_bf16 /proc/cpuinfo) $(date -u +%FT%TZ)"
    TIMEFORMAT='PROCESS_CPU user %U sys %S real %R'
    time taskset -c $core $C $extra -x "$(spec $k)" "$R/qaplib/wil50.dat"; } > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
wave() { one tri $1 0 & one box $1 1 & one tri $2 2 & one box $2 3 & wait; }
wave 10 13
wave d2 14
for k in 10 13 d2 14; do
  for arm in tri box; do
    f=$D/runs/${arm}_$k.log
    printf "%-4s %-4s %s | %s\n" $arm $k "$(grep '^quiet' $f | grep -o 'lb [0-9.]*')" "$(grep '^FINAL' $f)"
  done
done
