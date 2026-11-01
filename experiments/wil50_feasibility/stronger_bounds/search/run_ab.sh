#!/bin/bash
# The cut phase inside the search (bnb.C, QAP_DNN_CUT_ZONE), two A/Bs on
# 4 workers (bench/parallel.py), PRSM + degree 4, known optimum as incumbent:
#   nug30_{base,cut}   nug_deep.par -x 2:7:1:300:1; cut phase within 0.5% of
#                      the cutoff (31 units) with the 5x6 grid of the second
#                      matrix: the two depth-2 misses are 15 and 26 short
#   wil50n6_{base,cut} the subtree below the depth-3 near miss 20:12,40:28,17:25
#                      (node 6 of all32: 56 short at convergence, 82 at the
#                      give-up point), bnbdrv -F, DNN at depths 0-4 of the
#                      subtree; cut phase within 0.2% (98 units) with the 5x10
#                      grid. Base: the node gives up and its 47 children fathom;
#                      cut: the node itself should fathom.
# Logs: <name>.log; resumable (a log with "Total Time" is skipped).
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=4
run() {  # name instance ub extra [env...]
  f=$D/$1.log; grep -q "Total Time" "$f" 2>/dev/null && return
  shift; inst=$1; ub=$2; x=$3; shift 3
  (cd "$R" && env "$@" python3 bench/parallel.py qaplib/$inst.dat --jobs 4 --ub $ub \
     --params PARAM/nug_deep.par --extra "$x") > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
run wil50n6_cut  wil50 48816 "-F 20:12,40:28,17:25 -x 0:4:1:300:1" QAP_DNN_CUT_ZONE=0.002 QAP_DNN_CUT_GRID=5:10 QAP_DNN_CUT_SIDE=0
run wil50n6_base wil50 48816 "-F 20:12,40:28,17:25 -x 0:4:1:300:1" QAP_DNN_CUT_ZONE=0
run nug30_base   nug30 6124  "-x 2:7:1:300:1" QAP_DNN_CUT_ZONE=0
run nug30_cut    nug30 6124  "-x 2:7:1:300:1" QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_GRID=5:6 QAP_DNN_CUT_SIDE=1 QAP_DNN_CUT_DELTA=20
