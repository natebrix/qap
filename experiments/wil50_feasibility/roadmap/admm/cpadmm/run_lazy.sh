#!/bin/bash
# The lazy certificate inside the CP solve (DNN_CP_LAZYCERT=1, README.md "In the
# search"): node 13's separating CP solve at the production settings (run_sep.sh's
# s1sep100_prod_13 without the distance diagnostics), with and without the switch,
# ITERS (default 1500) iterations after the give-up point. Same bound at every check
# that certifies in both, the same final certified bound, less CPU per iteration.
#   bash run_lazy.sh CORE_OFF CORE_ON        logs: logs/cp_lazy{0,1}_13.log
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../../.." && pwd)
BIN=${BIN:-$R/experiments/dnn_cuts/dnncuts}
ITERS=${ITERS:-1500}
PROD_ENV="DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40 DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1"
cd "$R" || exit 1
mkdir -p "$D/logs"
one() { # lazy core
  local f=$D/logs/cp_lazy$1_13.log
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 $PROD_ENV DNN_CP_LAZYCERT=$1 taskset -c $2 "$BIN" -grid 5:10 -u 48816 -T 1 \
    -k 1000 -v -I 400 -x 21:28,8:41,28:19 -d 3 -cpsweeps 1 -cppasses 1 -cpsep 100:1000:0:0 -cpiters $ITERS -v \
    qaplib/wil50.dat > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
one 0 ${1:-0} &
one 1 ${2:-1} &
wait
