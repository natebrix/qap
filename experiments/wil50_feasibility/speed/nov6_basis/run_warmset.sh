#!/bin/bash
# The warm set of ../nov6_admm/run_warmset.sh (copied here so that its logs stay
# in this study's warm/): six wil50 depth-3 parents (given up at 300 iterations
# as in the search) and 24 of their depth-4 children, each child started from its
# parent's final Z, the conditioned Y (dnntest -K ... -Q b, the search's
# QAP_DNN_WARM_Y=2) and, with DNN_SHIP_BASIS=1, the parent's tracked basis.
#
# Flags: production (bundle D): degree 3, tight lambda (-L 2), PRSM (-g 1.0
# -G 1.0), refresh 40, give-up after 300 (-H 300), -i 1500 -u 48816, one
# CholQR pass, AMX, one BLAS thread, -v -v (the CPU profile of every solve:
# printing changes no number). One process per configuration, pinned to CPU.
#
# CONFIGS="name@env ..." (env: space-free assignments joined by '+', '-' for
# none); the production environment below is exported for every configuration
# (PROD=0: not). FAMILIES: names from the list below (default all six). EXTRA:
# more dnntest flags. BIN: the dnntest binary. OUT: the log directory (warm/).
# Resumable: a log with 5 RESULT lines is skipped. Logs: OUT/<name>_<family>.log
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
BIN=${BIN:-$R/latest/dnntest}
OUT=${OUT:-$D/warm}
CONFIGS=${CONFIGS:-"base@DNN_SHIP_BASIS=0 ship@DNN_SHIP_BASIS=1"}
FAMILIES=${FAMILIES:-"f9 f6 p10 p13 p2 p7"}
CPU=${CPU:-1}
# the production environment (pipeline_ab/run_pair.sh bestcore, the DNN part)
PRODENV="DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 DNN_ADMM_CROSS=2 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1 DNN_R32_NONCERT=1"
[ "${PROD:-1}" = 0 ] && PRODENV=""
declare -A PARENT KIDS
PARENT[f9]="26:27,8:38,34:28";   KIDS[f9]="37:14,37:37,37:21,37:26"
PARENT[f6]="20:12,40:28,17:25";  KIDS[f6]="43:18,43:7,43:22,43:19"
PARENT[p10]="39:27,30:48,25:49"; KIDS[p10]="7:5,26:5,3:5,47:5"
PARENT[p13]="21:28,8:41,28:19";  KIDS[p13]="3:2,44:2,15:2,18:2"
PARENT[p2]="18:43,11:29,46:44";  KIDS[p2]="8:21,9:21,30:21,37:21"
PARENT[p7]="12:7,33:6,32:11";    KIDS[p7]="5:34,8:34,42:34,0:34"
mkdir -p "$OUT"
for c in $CONFIGS; do
  name=${c%%@*}; envs=${c#*@}
  [ "$envs" = "-" ] && envs=""
  for fam in $FAMILIES; do
    f=$OUT/${name}_$fam.log
    [ "$(grep -c '^RESULT' "$f" 2>/dev/null)" = 5 ] && continue
    env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1 DNN_AMX=${AMX:-1} $PRODENV ${envs//+/ } \
      taskset -c "$CPU" "$BIN" -T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40 -H 300 -i 1500 -u 48816 -v -v \
      -x "${PARENT[$fam]}" -K "${KIDS[$fam]}" -Q b ${EXTRA} "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
  done
done
