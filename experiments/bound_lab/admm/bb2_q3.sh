#!/bin/bash
# B&B check of dnn.dscale on a solved sibling: tai25b (optimum 344,355,646 as the incumbent), tai50b's settings
# (bb_tai25b_base.par = experiments/estimates/tai50b/tai50b.par with dnn.depths = 1:23), without and with
# dnn.dscale = 0.5, on the depth-3 subtree of the optimal path (-F 0:3,1:14,2:9), 7-minute limit each,
# with the per-node DNN trace so that an unfinished run still gives per-node iterations.
cd "$(dirname "$0")"
OUT=$PWD/logs; PARD=$PWD
cd /home/user/wt-admm/latest
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for cfg in base ds0.5; do
  rm -f $OUT/bb_tai25b_F3_$cfg.dnn.txt $OUT/bb_tai25b_F3_$cfg.nodes.txt
  QAP_DNN_TRACE=$OUT/bb_tai25b_F3_$cfg.dnn.txt QAP_NODE_TRACE=$OUT/bb_tai25b_F3_$cfg.nodes.txt \
    timeout 420 ./bnbdrv -p $PARD/bb_tai25b_$cfg.par -u 344355646 -F 0:3,1:14,2:9 ../qaplib/tai25b.dat > $OUT/bb_tai25b_F3_$cfg.log 2>&1
  echo "exit $?" >> $OUT/bb_tai25b_F3_$cfg.log
done
