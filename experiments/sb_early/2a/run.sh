#!/bin/bash
# run.sh ARM TAG: wil50's depth-1 node 21:0 (first solve + cut phase) under the estimate's recipe with binary ARM
S=$S
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd /home/user/qap/latest && rm -f $S/m2a/$1.$2.trace && QAP_DNN_TRACE=$S/m2a/$1.$2.trace $S/bin/bnbdrv.$1 \
  -p ../infra/phase1/recipes/wil50_d2r10_g150xl2.par -s SYM/wil50t.sym -u 48818 -F 21:0 -d 1 -o $S/m2a/$1.$2.state ../qaplib/wil50t.dat > $S/m2a/$1.$2.log 2>&1
echo "$1 $2 exit $? $(date -u +%T)" >> $S/m2a/done.txt
