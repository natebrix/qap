#!/bin/bash
# A/B of the two settings ported from Jean-Pierre Goux's nug30 solver
# (QAP_DNN_EXT=20: the DNN bound kept up to 20 levels below a last-level
# miss; QAP_DNN_GU_ZONE=0.003: a give-up within 0.3% of the cutoff must be
# confirmed 100 iterations later), in the solver's PRSM + degree-4 mode,
# 4 workers, known optimum as the incumbent:
#   nug25  nug_deep.par -x 2:5:1:300:1   (give-up rule and warm starts)
#   nug28  nug_deep.par -x 2:5           (the setting where one depth-5 PRSM
#          miss sent 3.5 million nodes to QPB alone, eigensolver README sec. 8)
# Logs: <instance>_<base|goux>.log. Resumable: a log with "Total Time" is skipped.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=4
run() {  # name instance ub xflags [env...]
  f=$D/$1.log; grep -q "Total Time" "$f" 2>/dev/null && return
  shift; inst=$1; ub=$2; x=$3; shift 3
  (cd "$R" && env "$@" python3 bench/parallel.py qaplib/$inst.dat --jobs 4 --ub $ub \
     --params PARAM/nug_deep.par --extra "-x $x") > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
run nug25_base nug25 3744 2:5:1:300:1 QAP_DNN_EXT=0
run nug25_goux nug25 3744 2:5:1:300:1 QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003
run nug28_base nug28 5166 2:5 QAP_DNN_EXT=0
run nug28_goux nug28 5166 2:5 QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003
