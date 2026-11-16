#!/bin/bash
# allkeys/run.sh : every new key on at once, incumbent = optimum + 1 (+ 2 on nug: even objective values), so each
# run must find the optimum and prove it; a certified bound above a node's optimum on the optimal path would lose it.
H=$(cd "$(dirname "$0")" && pwd); cd /home/user/qap/latest
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for t in tai12b:39464925:1 tai15b:51765268:1 nug15:1150:2 nug16b:1240:2; do
  IFS=: read n opt d <<< "$t"
  QAP_DNN_POOL_CHECK=20 QAP_DNN_TRACE=$H/$n.dnn.txt QAP_NODE_TRACE=$H/$n.nodes.txt ./bnbdrv -p $H/$n.par -u $((opt + d)) ../qaplib/$n.dat > $H/$n.log 2> $H/$n.err
  echo "$n exit $? optimum $opt: $(grep 'New objective value' $H/$n.log | tail -n 1) | $(grep -E '^ *Nodes:' $H/$n.log) | $(grep -E '^ *(DNN bound|cut phase):' $H/$n.log | tr -s ' ' | tr '\n' ';') stderr $(wc -l < $H/$n.err) lines"
done
