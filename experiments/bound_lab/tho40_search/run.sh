#!/bin/bash
# One complete subtree, one arm, one solver process (cold from -F):
#   ./run.sh NAME ARM SPEC      -> runs/NAME_ARM/{log.txt,nodes.txt,dnn.txt,cmd.txt}
# ARM names par/ARM.par. Instance, symmetry file and incumbent as experiments/estimates/tho40/settings.sh.
set -u
H=$(cd "$(dirname "$0")" && pwd); R=/home/user/qap
NAME=$1; ARM=$2; SPEC=$3
D=$H/runs/${NAME}_$ARM; mkdir -p "$D"
[ -f "$D/done" ] && { echo "$D already done"; exit 0; }
rm -f "$D/nodes.txt" "$D/dnn.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt
cd $R/latest
CMD="./bnbdrv -p $H/par/$ARM.par -s SYM/tho40t.sym -u 240516 -F $SPEC ../qaplib/tho40t.dat"
{ echo "$CMD"; echo "QAP_NODE_TRACE=$QAP_NODE_TRACE QAP_DNN_TRACE=$QAP_DNN_TRACE OPENBLAS_NUM_THREADS=1"; date -u +%FT%TZ; uptime; } > "$D/cmd.txt"
$CMD > "$D/log.txt" 2> "$D/err.txt"
echo "exit $?" >> "$D/cmd.txt"; date -u +%FT%TZ >> "$D/cmd.txt"; uptime >> "$D/cmd.txt"
touch "$D/done"
