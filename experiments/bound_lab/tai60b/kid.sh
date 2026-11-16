#!/bin/bash
# kid.sh TAG PAR INSTATE DEPTH F:L [UB]: one node, warm from its parent's saved state (after regions/size50/kid.sh).
# INSTATE (a parent's -o state, its children set aside by -d) is cut to the depth-DEPTH child whose permutation
# has entry F equal to L (tho40_search/cutstate.py); that node is bounded under par/PAR.par (first solve with the
# inherited pool, then the phase) and, if it does not fathom, its children are written to runs/TAG/out.state
# (-d DEPTH -o) for the next level. UB: another incumbent (validity runs at best known + 1).
D=$(cd $(dirname $0); pwd); I=${INST:-tai60b}
TAG=$1; PAR=$2; IN=$3; DEPTH=$4; FL=$5; UB=${6:-}
mkdir -p $D/runs/$TAG; rm -f $D/runs/$TAG/{nodes,dnn,live}.txt
python3 -I $D/../tho40_search/cutstate.py $IN $D/runs/$TAG/in.state $DEPTH $FL || { echo "cutstate failed $TAG"; exit 1; }
if [ -n "$UB" ]; then sed -i "s/^incumbent [0-9.e+]* /incumbent $UB /" $D/runs/$TAG/in.state; fi
cd /home/user/wt-regions/latest
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 QAP_NODE_TRACE=$D/runs/$TAG/nodes.txt QAP_DNN_TRACE=$D/runs/$TAG/dnn.txt QAP_DNN_TRACE_HIST=1 QAP_CP_LIVE=$D/runs/$TAG/live.txt
CMD="./bnbdrv -p $D/par/$PAR.par --allow-settings-change -r $D/runs/$TAG/in.state -d $DEPTH -o $D/runs/$TAG/out.state ../qaplib/$I.dat"
{ echo "$CMD"; date -u +%FT%TZ; uptime; } > $D/runs/$TAG/cmd.txt
$CMD > $D/runs/$TAG/log.txt 2> $D/runs/$TAG/err.txt
{ echo "exit $?"; date -u +%FT%TZ; uptime; } >> $D/runs/$TAG/cmd.txt
rm -f $D/runs/$TAG/in.state
grep -q "^node " $D/runs/$TAG/out.state 2>/dev/null || rm -f $D/runs/$TAG/out.state   # fathomed: nothing to keep
cut -c1-170 $D/runs/$TAG/nodes.txt
