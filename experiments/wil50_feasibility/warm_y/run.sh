#!/bin/bash
# Full-state child warm start (roadmap/admm/splitting.md idea 1): does a child
# start faster from its parent's final Y as well as its Z? The 16 depth-4
# children of stronger_bounds/depth4/nodes.txt, each bounded after its depth-3
# parent (a clear miss, given up at 300 iterations as in the search), from
# three starts (dnntest -Q):
#   z  the parent's Z only (today's search)
#   a  Z and the parent's Y restricted to the child's indices
#   b  Z and the parent's Y conditioned on the child's new assignment
# Production settings (bundle D, amx/README.md): DNN_AMX=1, one-pass CholQR,
# degree 3 with tight lambda, refresh 40, PRSM (1, 1), give-up after 300,
# at most 1,500 iterations. 4 at a time; resumable (a log with three RESULT
# lines for the child is skipped).
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../.." && pwd)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_AMX=1 DNN_CHOLQR_PASSES=1
one() {
  k=$1; spec=$2
  f=$D/runs/child_$k.log
  [ "$(grep -c '^RESULT node=[123] ' "$f" 2>/dev/null)" = 3 ] && exit 0
  parent=${spec%,*}; child=${spec##*,}
  "$R/latest/dnntest" -T 1 -d 3 -L 2 -r 40 -g 1.0 -G 1.0 -H 300 -i 1500 -u 48816 -v \
    -x "$parent" -K "$child" -Q z,a,b "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
export -f one; export D R
awk '{print $1, $3}' "$R/experiments/wil50_feasibility/stronger_bounds/depth4/nodes.txt" | xargs -P 4 -L 1 bash -c 'one "$@"' _
