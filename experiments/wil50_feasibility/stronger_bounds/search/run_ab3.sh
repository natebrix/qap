#!/bin/bash
# Second attempt at the dominant-class end-to-end test. run_ab2.sh was stopped
# after 8 h: with the DNN bound limited to four levels below the node
# (-x 0:4), every miss at the last DNN level handed its subtree to the QPB
# bound alone, and at 43 free facilities those subtrees ran past 10^7 QPB nodes
# per worker (partial_ab2/). Same subtree, the depth-3 clear miss node 10 of
# all32 (39:27,30:48,25:49; 158 short at convergence), same settings, plus
# QAP_DNN_EXT=20 (the DNN keeps running below a last-level miss, Goux D8) and
# a parallel.py checkpoint every 2 hours (a checkpoint idles the early
# finishers until the slowest node ends, up to a cut phase: 10 minutes was
# costing a quarter of the machine) and a progress line per finished task.
# Cut arm: zone 0.003, quiet round of 1,500 from the best round's multipliers,
# stall rule 3. Base: no cut phase. 4 workers, PRSM + degree 4. Resumable:
# a finished arm has a "Total Time" line; an unfinished one resumes from its
# checkpoint.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=4 QAP_DNN_EXT=20
run() {  # name instance ub extra [env...]
  name=$1; f=$D/$name.log; grep -q "Total Time" "$f" 2>/dev/null && return
  shift; inst=$1; ub=$2; x=$3; shift 3
  resume=(); [ -s "$D/$name.ckpt" ] && resume=(--resume "$D/$name.ckpt")
  (cd "$R" && env "$@" python3 bench/parallel.py qaplib/$inst.dat --jobs 4 --ub $ub \
     --params PARAM/nug_deep.par --extra "$x" --checkpoint "$D/$name.ckpt" --checkpoint-every 7200 "${resume[@]}") \
     >> "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
run wil50n10ext_cut  wil50 48816 "-F 39:27,30:48,25:49 -x 0:4:1:300:1" QAP_DNN_CUT_ZONE=0.003 QAP_DNN_CUT_GRID=5:10 QAP_DNN_CUT_SIDE=0 QAP_DNN_CUT_FREEZE=1500
run wil50n10ext_base wil50 48816 "-F 39:27,30:48,25:49 -x 0:4:1:300:1" QAP_DNN_CUT_ZONE=0
