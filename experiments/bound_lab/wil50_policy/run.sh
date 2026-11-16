#!/bin/bash
# One complete subtree, one arm, one solver process (-T 1), cold from -F:
#   ./run.sh INST NAME ARM SPEC [EXTRA]   -> runs/NAME_ARM/{cmd.txt,log.txt,err.txt,nodes.txt,dnn.txt}
# INST wil50 | sko49 (instance, symmetry file, incumbent as the instance's settings.sh);
# ARM names par/INST_ARM.par. NAME w*: warm, resumes runs/NAME_ARM/in.state (made by warm.sh).
# The solver is latest/bnbdrv with tho40_search/stop/handover.patch applied (unset key = unpatched).
set -u
H=$(cd "$(dirname "$0")" && pwd); R=$(cd "$H/../../.." && pwd)
INST=$1; NAME=$2; ARM=$3; SPEC=$4; EXTRA=${5:-}
case $INST in wil50) UB=48816;; sko49) UB=23386;; *) echo "INST?"; exit 2;; esac
D=$H/runs/${NAME}_$ARM; mkdir -p "$D"
[ -f "$D/done" ] && { echo "$D already done"; exit 0; }
rm -f "$D/nodes.txt" "$D/dnn.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt QAP_DNN_TRACE_HIST=1 QAP_DNN_POOL_CHECK=20
START="-F $SPEC"
case $NAME in w*) START="-r $D/in.state --allow-settings-change";; esac
cd "$R/latest"
CMD="./bnbdrv -p $H/par/${INST}_$ARM.par -s SYM/${INST}t.sym -u $UB $START $EXTRA $R/qaplib/${INST}t.dat"
{ echo "$CMD"; echo "QAP_NODE_TRACE QAP_DNN_TRACE QAP_DNN_TRACE_HIST=1 QAP_DNN_POOL_CHECK=20 OPENBLAS_NUM_THREADS=1"; date -u +%FT%TZ; uptime; } > "$D/cmd.txt"
$CMD > "$D/log.txt" 2> "$D/err.txt"
echo "exit $?" >> "$D/cmd.txt"; date -u +%FT%TZ >> "$D/cmd.txt"; uptime >> "$D/cmd.txt"
touch "$D/done"
