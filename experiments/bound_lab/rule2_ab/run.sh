#!/bin/bash
# run.sh NAME ARM FIXES : one tho40 subtree under one par, traced
D=/home/user/qap/experiments/bound_lab/rule2_ab; R=$D/runs/$1_$2; mkdir -p $R
cd /home/user/qap/latest
export QAP_NODE_TRACE=$R/nodes.txt QAP_DNN_TRACE=$R/dnn.txt OPENBLAS_NUM_THREADS=1
echo "./bnbdrv -p $D/$2.par -s SYM/tho40t.sym -u 240516 -F $3 ../qaplib/tho40t.dat" > $R/cmd.txt
./bnbdrv -p $D/$2.par -s SYM/tho40t.sym -u 240516 -F $3 ../qaplib/tho40t.dat > $R/log.txt 2>&1
echo done > $R/done
