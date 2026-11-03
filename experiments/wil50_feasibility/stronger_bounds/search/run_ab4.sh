#!/bin/bash
# Branching at the DNN depths (roadmap/tree_size.md item 2): the node-10 subtree
# cut arm of run_ab3 (nug_deep.par: rules 4/7/6 by gap, which sends the DNN
# depths to rule 2; 365 nodes, 76,755 CPU-s, wil50n10ext_cut.log) against the
# same arm branching on the DNN marginal lines at every node through depth 6:
#   rule9  PARAM/nug_deep_xbar9_all.par  (lookahead over rule 8's two marginal
#          lines, the children's QPB bounds choose; the hedge)
#   rule8  PARAM/nug_deep_xbar_all.par   (the line with the largest sum of
#          squared marginals)
# Everything else as run_ab3's cut arm: -x 0:4:1:300:1, QAP_DNN_EXT=20, zone
# 0.003, grid 5:10 side 0, quiet 1,500, PRSM + degree 4, 4 workers, checkpoint
# every 2 h. Resumable (Total Time line = done; else from the .ckpt).
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=4 QAP_DNN_EXT=20
export QAP_DNN_CUT_ZONE=0.003 QAP_DNN_CUT_GRID=5:10 QAP_DNN_CUT_SIDE=0 QAP_DNN_CUT_FREEZE=1500
run() {  # name params
  name=$1; par=$2; f=$D/$name.log; grep -q "Total Time" "$f" 2>/dev/null && return
  resume=(); [ -s "$D/$name.ckpt" ] && resume=(--resume "$D/$name.ckpt")
  (cd "$R" && python3 bench/parallel.py qaplib/wil50.dat --jobs 4 --ub 48816 \
     --params $par --extra "-F 39:27,30:48,25:49 -x 0:4:1:300:1" --checkpoint "$D/$name.ckpt" --checkpoint-every 7200 "${resume[@]}") \
     >> "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
run wil50n10ext_cut_rule9 PARAM/nug_deep_xbar9_all.par
run wil50n10ext_cut_rule8 PARAM/nug_deep_xbar_all.par
