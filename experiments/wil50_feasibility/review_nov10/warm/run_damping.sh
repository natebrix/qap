#!/bin/bash
# Warm children of two depth-3 parents (families f9 and f6 of
# ../../speed/nov6_admm/run_warmset.sh), a certificate at every iteration, at
# the full run's settings. Extra environment after the name, e.g.
#   run_damping.sh 0 base
#   run_damping.sh 1 d2 DNN_ADMM_WARM_DAMP=1e-9:2
# Logs: NAME_f9.log, NAME_f6.log here; read them with rd.py f9 f6.
W=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$W/../../../.." && pwd)
CPU=$1; NAME=$2; shift 2
declare -A PARENT KIDS
PARENT[f9]="26:27,8:38,34:28"; KIDS[f9]="37:14,37:37,37:21,37:26"
PARENT[f6]="20:12,40:28,17:25"; KIDS[f6]="43:18,43:7,43:22,43:19"
for fam in ${FAMS:-f9 f6}; do
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "$@" taskset -c "$CPU" "$R/latest/dnntest" \
    -p "$R/experiments/wil50_feasibility/fullrun/fullrun.par" -T 1 -i 1500 -u 48816 -v ${KFLAG--k 1} \
    -x "${PARENT[$fam]}" -K "${KIDS[$fam]}" -Q b "$R/qaplib/wil50.dat" > "$W/${NAME}_${fam}.log" 2>&1
done
