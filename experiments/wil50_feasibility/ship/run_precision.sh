#!/bin/bash
# Which precision the shipped parent state needs on wil50 (QAP_DNN_SHIP):
# node 10 bounded for 300 iterations as in the search, then 8 of its depth-4
# children started from its Z and conditioned Y (-Q b, QAP_DNN_WARM_Y=2),
# with the parent state exact, rounded to fp32 (-Wr 32) or to bf16 (-Wr 16),
# as dnntest's -Wr does (the convergence study's test S, on more children).
# Three runs at once on cores 0-2; logs in logs/.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../.." && pwd)
DT=$R/latest/dnntest; W=$R/qaplib/wil50.dat
mkdir -p "$D/logs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_AMX=1 DNN_CHOLQR_PASSES=1
KIDS=7:45,47:45,37:45,17:45,27:45,46:45,6:45,36:45
core=0
for wr in exact 32 16; do
  extra=(); [ "$wr" != exact ] && extra=(-Wr "$wr")
  taskset -c $core "$DT" -T 1 -d 3 -L 2 -r 40 -g 1.0 -G 1.0 -H 300 -i 1500 -u 48816 \
    -x 39:27,30:48,25:49 -J 300 "${extra[@]}" -K "$KIDS" -Q b "$W" > "$D/logs/n10_kids_$wr.log" 2>&1 &
  core=$((core + 1))
done
wait
for wr in exact 32 16; do echo "== $wr"; grep -E "^CHILD-SPEC|^RESULT" "$D/logs/n10_kids_$wr.log" | paste - - | awk '{print $2, $6, $8, $5}'; done
