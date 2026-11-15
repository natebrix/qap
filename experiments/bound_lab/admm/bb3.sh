#!/bin/bash
# B&B check of dnn.dscale on tai30b (optimum 637,117,113 as the incumbent; root gap 6%, like tai50b's):
# the depth-3 subtree of the optimal path (-F 0:3,1:7,2:10), tai50b's settings with dnn.depths = 1:28,
# without and with dnn.dscale = 0.5, 9-minute limit each, per-node traces.
cd "$(dirname "$0")"


OUT=$PWD/logs; PARD=$PWD
cd /home/user/wt-admm/latest
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for cfg in base ds0.5; do
  rm -f $OUT/bb_tai30b_F3_$cfg.dnn.txt $OUT/bb_tai30b_F3_$cfg.nodes.txt
  QAP_DNN_TRACE=$OUT/bb_tai30b_F3_$cfg.dnn.txt QAP_NODE_TRACE=$OUT/bb_tai30b_F3_$cfg.nodes.txt \
    timeout 540 ./bnbdrv -p $PARD/bb_tai30b_$cfg.par -u 637117113 -F 0:3,1:7,2:10 ../qaplib/tai30b.dat > $OUT/bb_tai30b_F3_$cfg.log 2>&1
  echo "exit $?" >> $OUT/bb_tai30b_F3_$cfg.log
done
