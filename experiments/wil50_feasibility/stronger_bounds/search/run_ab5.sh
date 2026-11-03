#!/bin/bash
# Cut-phase trims (roadmap/README.md item 3, wil50_feasibility/trims/): the
# rule-9 arm of run_ab4 (node-10 subtree, PARAM/nug_deep_xbar9_all.par, cut
# phase zone 0.003, quiet 1,500, stall 3; 322 nodes, 54,621 CPU-s,
# wil50n10ext_cut_rule9.log) against the same arm with
#   QAP_DNN_CUT_PGATE=4:20  the primal gate: a phase stops after round 4 when
#                           its best primal is more than 20 below the cutoff
#                           (replay on 88 logged phases: 44 of 63 misses
#                           stopped, no false stop; trims/replay_pgate.py)
#   QAP_DNN_CUT_ZONE=0.0035 the wider zone the tree memo's model prefers once
#                           the gate makes a failed phase cheap
# The give-up deferral (QAP_DNN_GU_ZONE) stays off: the replay on the depth-4
# children (trims/replay_giveup.py) shows the extrapolation rule already lets
# the slow uncut fathoms finish; the deferral only adds cost.
# Everything else as run_ab4: -x 0:4:1:300:1, QAP_DNN_EXT=20, grid 5:10 side 0,
# PRSM + degree 4, 4 workers, checkpoint every 2 h. Resumable.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=4 QAP_DNN_EXT=20
export QAP_DNN_CUT_ZONE=0.0035 QAP_DNN_CUT_GRID=5:10 QAP_DNN_CUT_SIDE=0 QAP_DNN_CUT_FREEZE=1500
export QAP_DNN_CUT_PGATE=4:20
run() {  # name params
  name=$1; par=$2; f=$D/$name.log; grep -q "Total Time" "$f" 2>/dev/null && return
  resume=(); [ -s "$D/$name.ckpt" ] && resume=(--resume "$D/$name.ckpt")
  (cd "$R" && python3 bench/parallel.py qaplib/wil50.dat --jobs 4 --ub 48816 \
     --params $par --extra "-F 39:27,30:48,25:49 -x 0:4:1:300:1" --checkpoint "$D/$name.ckpt" --checkpoint-every 7200 "${resume[@]}") \
     >> "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
run wil50n10ext_cut_rule9_trim PARAM/nug_deep_xbar9_all.par
