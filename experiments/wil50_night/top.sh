#!/bin/bash
# wil50_night/top.sh ARM [NODE]: the line an arm picks at the root (no NODE) or at a depth-1 node, cold
# (-d stops the run once the node has branched). Output in top/ARM_NODE.{log,nodes.txt}.
D=$(cd "$(dirname "$0")"; pwd); BUILD=${BUILD:-/home/user/wt-rule10}
arm=$1; node=${2:-}; tag=${arm}_${node:-root}; tag=${tag//:/-}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$BUILD/latest" && QAP_NODE_TRACE="$D/top/$tag.nodes.txt" ./bnbdrv -p "$D/par/$arm.par" -s SYM/wil50t.sym -u 48818 \
  ${node:+-F "$node"} -d $([ -n "$node" ] && echo 1 || echo 0) -o "$D/top/$tag.state" ../qaplib/wil50t.dat > "$D/top/$tag.log" 2> "$D/top/$tag.err"
echo "exit $? $(date -u +%T)" >> "$D/top/$tag.log"
