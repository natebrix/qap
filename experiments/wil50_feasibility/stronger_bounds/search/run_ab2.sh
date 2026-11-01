#!/bin/bash
# The dominant class end to end: the subtree below the depth-3 clear miss
# node 10 of all32 (39:27,30:48,25:49; 158 short at convergence). Base: the
# node gives up, 47 depth-4 children of which about 44% fathom uncut and the
# rest spawn 46 depth-5 children each. Cut: zone 0.003 (146 units at the
# give-up point), the quiet round from the best round's multipliers with 1,500
# iterations, stall rule 3. 4 workers, PRSM + degree 4. Resumable.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=4
run() {  # name instance ub extra [env...]
  f=$D/$1.log; grep -q "Total Time" "$f" 2>/dev/null && return
  shift; inst=$1; ub=$2; x=$3; shift 3
  (cd "$R" && env "$@" python3 bench/parallel.py qaplib/$inst.dat --jobs 4 --ub $ub \
     --params PARAM/nug_deep.par --extra "$x") > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
run wil50n10_cut  wil50 48816 "-F 39:27,30:48,25:49 -x 0:4:1:300:1" QAP_DNN_CUT_ZONE=0.003 QAP_DNN_CUT_GRID=5:10 QAP_DNN_CUT_SIDE=0 QAP_DNN_CUT_FREEZE=1500
run wil50n10_base wil50 48816 "-F 39:27,30:48,25:49 -x 0:4:1:300:1" QAP_DNN_CUT_ZONE=0
