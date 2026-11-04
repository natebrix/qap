#!/bin/bash
# E4 (question 3): free child bounds at a give-up. Node 10 (39:27,30:48,25:49)
# with the production flags to its give-up at 300 iterations, M = L + scale Zc
# of the check at iteration 300 dumped (DNN_CERTX bit 16; the dump, 39 MB,
# goes to the scratch directory given), -C 2 for the reduced-cost child bounds
# of the standard certificate; then scripts/child_gl.py on the branching line
# of the search (column 45, the rule-9 line in pool_search/trace.txt).
# Usage: run_e4.sh <core> <dumpdir>
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
core=$1; dd=$2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export DNN_AMX=1 DNN_CHOLQR_PASSES=1
f=$D/logs/e4_n10.log
if ! grep -q '^RESULT' "$f" 2>/dev/null; then
  DNN_CERTX=17 DNN_CERTX_EVERY=30 DNN_CERTX_WINDOW=-1e9 DNN_CERTX_DUMP=$dd/n10 \
    taskset -c "$core" "$R/latest/dnntest" -T 1 -g 1.0 -G 1.0 -d 3 -L 2 -r 40 -H 300 -u 48816 -C 2 \
    -x "39:27,30:48,25:49" "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
fi
# the local column of location 45 at node 10 (locations 27, 48, 49 are fixed)
taskset -c "$core" python3 "$D/scripts/child_gl.py" "$dd/n10" line:col:44 > "$D/logs/e4_n10_childgl_col45.log" 2>&1
