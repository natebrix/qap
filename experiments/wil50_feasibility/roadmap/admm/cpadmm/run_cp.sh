#!/bin/bash
# CP-ADMM on node 13 (21:28,8:41,28:19): the round-12 pool of run_ref.sh
# (logs/pool12_13.txt, 13,000 regional triangles, frozen) inside the Y-step, from
# the give-up point (the uncut round 0 of the rounds, 400 iterations), certified at
# every check (every 10 iterations) with the Dykstra multipliers. The dual test's
# settings (dnncuts defaults: -d 4, refresh 20, gamma 1/1, beta 0.3) unless the
# arm says prod (the production DNN settings and bundle, see below).
#   bash run_cp.sh ARM [ARM...]     ARM = name:core; one job per arm, in parallel
# Arms (logs/cp_<name>_13.log; resumable, FINAL line = done):
#   s1p1   1 Dykstra sweep per iteration, 1 half-space pass per sweep (the memo's design)
#   s1p5   1 sweep of 5 passes
#   s3p3   3 sweeps of 3 passes
#   plain  -cpplain: the same solve without the cuts (uncut continuation; the
#          plain PRSM iteration's cost under the same load)
#   *_prod the same with the production settings: DNN_AMX=1 DNN_CHOLQR_PASSES=1,
#          -d 3, tight lambda 2, refresh 40, gamma 1/1, and DNN_R32_NONCERT=1
#          DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 DNN_CORE_LAZYCERT=1
#          DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1 (LAZYCERT is
#          ignored by the CP solve, which certifies every check)
# ITERS (default 2000) iterations after the give-up point; DIAG (default 200:30):
# every 200 iterations the distance of the Y-step's output to a 30-sweep projection;
# QUIET (default 0): then a quiet round of that many iterations with mu_D frozen.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../../.." && pwd)
BIN=${BIN:-$R/experiments/dnn_cuts/dnncuts}
POOL=${POOL:-$D/logs/pool12_13.txt}
ITERS=${ITERS:-2000}
DIAG=${DIAG:-200:30}
QUIET=${QUIET:-0}
TAG=${TAG:-}
BASE="-grid 5:10 -u 48816 -T 1 -v -I 400 -x 21:28,8:41,28:19"
PROD_ENV="DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40 DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1"
cd "$R" || exit 1
mkdir -p "$D/logs"
one() { # name core
  local name=$1 core=$2 envs="" flags="" b=${1%_prod}
  [ "$b" != "$1" ] && { envs=$PROD_ENV; flags="-d 3"; }
  case $b in
    s1p1) flags="$flags -cpadmm $POOL -cpsweeps 1 -cppasses 1 -cpdiag $DIAG" ;;
    s1p5) flags="$flags -cpadmm $POOL -cpsweeps 1 -cppasses 5 -cpdiag $DIAG" ;;
    s3p3) flags="$flags -cpadmm $POOL -cpsweeps 3 -cppasses 3 -cpdiag $DIAG" ;;
    s1p20) flags="$flags -cpadmm $POOL -cpsweeps 1 -cppasses 20 -cpdiag $DIAG" ;;
    plain) flags="$flags -cpplain" ;;
    *) echo "unknown arm $name"; return 1 ;;
  esac
  local f=$D/logs/cp_${name}${TAG}_13.log
  grep -q '^FINAL' "$f" 2>/dev/null && return 0
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 $envs taskset -c $core "$BIN" $BASE $flags -cpiters $ITERS -cpquiet $QUIET -v qaplib/wil50.dat > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
for a in "$@"; do one ${a%:*} ${a#*:} & done
wait
