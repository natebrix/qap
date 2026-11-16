#!/bin/bash
# root.sh NAME [INST UB]: the root of INST (default tai60b, incumbent 608215054 = best known) in the search under
# par/NAME.par (DNN at depth 0, then the cut phase), its children set aside (-d 0) and written with the root's
# shipped DNN state and pool to runs/NAME/root.state. runs/NAME/live.txt: the phase's running best at each check.
D=$(cd $(dirname $0); pwd); N=$1; I=${2:-tai60b}; UB=${3:-608215054}; mkdir -p $D/runs/$N; cd /home/user/wt-regions/latest
rm -f $D/runs/$N/{nodes,dnn,live}.txt
export OPENBLAS_NUM_THREADS=${THREADS:-1} OMP_NUM_THREADS=1 QAP_NODE_TRACE=$D/runs/$N/nodes.txt QAP_DNN_TRACE=$D/runs/$N/dnn.txt QAP_DNN_TRACE_HIST=1 QAP_CP_LIVE=$D/runs/$N/live.txt
CMD="./bnbdrv -p $D/par/$N.par -u $UB -d 0 -o $D/runs/$N/root.state ../qaplib/$I.dat"
{ echo "$CMD"; date -u +%FT%TZ; uptime; } > $D/runs/$N/cmd.txt
$D/rss.sh $D/runs/$N/log.txt $CMD
{ echo "exit $?"; tail -1 $D/runs/$N/log.txt; date -u +%FT%TZ; uptime; } >> $D/runs/$N/cmd.txt
