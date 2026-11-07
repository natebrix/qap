#!/bin/bash
# CP-ADMM with cuts separated inside the solve (README.md, "Separation inside the
# solve"). From the give-up point (round 0: 400 uncut iterations, then its own
# separation of 1,000 triangles, the start pool), one CP solve that separates up to
# 1,000 new violated regional triangles on its projected Y every k iterations
# (-cpsep k:1000[:CAP[:DROP]]). The dual test's settings unless the arm ends in _prod
# (run_cp.sh's production settings).
#   bash run_sep.sh ARM:CORE [ARM:CORE ...]       one job per arm, in parallel
# ARM = <sN>sep<k>[_prod]_<node>: s1 = 1 sweep of 1 pass, s3 = 3 sweeps of 3 passes;
# k the cadence; node 13, d2 or 14. Examples: s1sep300_13, s1sep100_d2, s3sep100_13,
# s1sep100_prod_13. Log: logs/cp_<ARM>.log (resumable; FINAL line = done).
# ITERS (default 3000) iterations after the give-up point; CAP, DROP (default 0: off).
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../../.." && pwd)
BIN=${BIN:-$R/experiments/dnn_cuts/dnncuts}
ITERS=${ITERS:-3000}
CAP=${CAP:-0}
DROP=${DROP:-0}
PROD_ENV="DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40 DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1"
cd "$R" || exit 1
mkdir -p "$D/logs"
spec() { case $1 in 13) echo 21:28,8:41,28:19;; d2) echo 39:27,30:48,25:49,26:5;; 14) echo 47:23,10:44,40:34;; esac; }
one() { # arm core
  local arm=$1 core=$2 node=${1##*_} envs="" flags=""
  local body=${arm%_*} sw k
  case $body in *_prod) envs=$PROD_ENV; flags="-d 3"; body=${body%_prod};; esac
  sw=${body%%sep*}; k=${body#*sep}
  case $sw in s1) flags="$flags -cpsweeps 1 -cppasses 1";; s3) flags="$flags -cpsweeps 3 -cppasses 3";; *) echo "bad arm $arm"; return 1;; esac
  local x; x=$(spec $node); [ -z "$x" ] && { echo "bad node in $arm"; return 1; }
  local f=$D/logs/cp_${arm}.log
  grep -q '^FINAL' "$f" 2>/dev/null && return 0
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 $envs taskset -c $core "$BIN" -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 \
    -x $x $flags -cpsep $k:1000:$CAP:$DROP -cpdiag 500:30 -cpiters $ITERS -v qaplib/wil50.dat > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
for a in "$@"; do one ${a%:*} ${a#*:} & done
wait
