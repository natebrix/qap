#!/bin/bash
# root.sh NAME: the tai50b root in the search under par/NAME.par (DNN at depth 0, cut phase), its children
# set aside (-d 0) and written with the root's shipped DNN state and pool to runs/NAME/root.state
D=/home/user/qap/experiments/bound_lab/regions/size50; N=$1; mkdir -p $D/runs/$N; cd /home/user/wt-regions/latest
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 QAP_NODE_TRACE=$D/runs/$N/nodes.txt QAP_DNN_TRACE=$D/runs/$N/dnn.txt QAP_DNN_TRACE_HIST=1
CMD="./bnbdrv -p $D/par/$N.par -u 458821517 -d 0 -o $D/runs/$N/root.state ../qaplib/tai50b.dat"
{ echo "$CMD"; date -u +%FT%TZ; uptime; } > $D/runs/$N/cmd.txt
$CMD > $D/runs/$N/log.txt 2> $D/runs/$N/err.txt
{ echo "exit $?"; date -u +%FT%TZ; uptime; } >> $D/runs/$N/cmd.txt
