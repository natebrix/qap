#!/bin/zsh
# The cut phase in the search on the four wil50 depth-3 near misses (all32 ids
# 6, 9, 15, 16), laptop: base (the node gives up, its 47 children fathom),
# cut phase as the cloud ran it (cap at the cutoff), and cut phase with
# QAP_DNN_CUT_NOCAP=1 QAP_DNN_CUT_FREEZE=600. Same solver settings as
# ../../search/run_ab.sh: PRSM + degree 4, nug_deep.par, -x 0:4:1:300:1.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../../.." && pwd)
export VECLIB_MAXIMUM_THREADS=1 QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=4
typeset -A N; N[6]=20:12,40:28,17:25; N[9]=26:27,8:38,34:28; N[15]=49:20,5:16,23:28; N[16]=30:37,27:27,49:16
run() { # name jobs node env...
  f=$D/$1.log; grep -q "Total Time" "$f" 2>/dev/null && return
  j=$2; node=$3; shift 3
  (cd $R/latest && env "$@" python3 ../bench/parallel.py ../qaplib/wil50.dat --jobs $j --ub 48816 \
     --params PARAM/nug_deep.par --extra "-F $node -x 0:4:1:300:1") > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
CUT=(QAP_DNN_CUT_ZONE=0.002 QAP_DNN_CUT_GRID=5:10 QAP_DNN_CUT_SIDE=0)
( for k in 6 9 15 16; do run n${k}_cut   1 $N[$k] $CUT; run n${k}_nocap 1 $N[$k] $CUT QAP_DNN_CUT_NOCAP=1 QAP_DNN_CUT_FREEZE=600; done ) &
( for k in 6 9 15 16; do run n${k}_base  4 $N[$k] QAP_DNN_CUT_ZONE=0; done ) &
wait
