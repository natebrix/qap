#!/bin/bash
# as bb3.sh (tai30b, the depth-3 subtree of the optimal path), three more arms: dnn.dscale = 0.25, and
# dnn.dscale = 0.5 for the first 200 or 100 iterations of each call (dnn.dscale_until), the plain iteration after
cd "$(dirname "$0")"
while pgrep -f "q8.sh" > /dev/null; do sleep 5; done
OUT=$PWD/logs; PARD=$PWD
cd /home/user/wt-admm/latest
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for cfg in ds0.5u200 ds0.25 ds0.5u100; do
  rm -f $OUT/bb_tai30b_F3_$cfg.dnn.txt $OUT/bb_tai30b_F3_$cfg.nodes.txt
  QAP_DNN_TRACE=$OUT/bb_tai30b_F3_$cfg.dnn.txt QAP_NODE_TRACE=$OUT/bb_tai30b_F3_$cfg.nodes.txt \
    timeout 540 ./bnbdrv -p $PARD/bb_tai30b_$cfg.par -u 637117113 -F 0:3,1:7,2:10 ../qaplib/tai30b.dat > $OUT/bb_tai30b_F3_$cfg.log 2>&1
  echo "exit $?" >> $OUT/bb_tai30b_F3_$cfg.log
done
