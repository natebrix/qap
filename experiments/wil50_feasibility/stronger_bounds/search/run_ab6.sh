#!/bin/bash
# The stacked A/B (Nov 3): the run_ab4 rule-9 arm (node-10 subtree, cut phase
# zone 0.003, quiet 1,500, stall 3, EXT 20, rule 9 at every node through
# depth 6; 322 nodes, 54,621 CPU-s, 4 h 22 min, wil50n10ext_cut_rule9.log)
# with every per-node lever measured since, all opt-in switches:
#   bundle D (amx/README.md): DNN_AMX=1 DNN_CHOLQR_PASSES=1, degree 3 with
#     tight lambda (QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2), refresh every 40
#     (QAP_DNN_REFRESH=40); 1.71x per node on the 32-node sample
#   children start from the parent's Y as well as Z, conditioned on the new
#     assignment (QAP_DNN_WARM_Y=2, warm_y/README.md); 1.17x on 16 children
# Everything else as run_ab4: PRSM (1, 1), -x 0:4:1:300:1, 4 workers,
# checkpoint every 2 h. Resumable (Total Time line = done; else the .ckpt).
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_EXT=20
export QAP_DNN_CUT_ZONE=0.003 QAP_DNN_CUT_GRID=5:10 QAP_DNN_CUT_SIDE=0 QAP_DNN_CUT_FREEZE=1500
export DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40
export QAP_DNN_WARM_Y=2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
run() {  # name params
  name=$1; par=$2; f=$D/$name.log; grep -q "Total Time" "$f" 2>/dev/null && return
  resume=(); [ -s "$D/$name.ckpt" ] && resume=(--resume "$D/$name.ckpt")
  (cd "$R" && python3 bench/parallel.py qaplib/wil50.dat --jobs 4 --ub 48816 \
     --params $par --extra "-F 39:27,30:48,25:49 -x 0:4:1:300:1" --checkpoint "$D/$name.ckpt" --checkpoint-every 7200 "${resume[@]}") \
     >> "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
run wil50n10ext_cut_rule9_stacked PARAM/nug_deep_xbar9_all.par
