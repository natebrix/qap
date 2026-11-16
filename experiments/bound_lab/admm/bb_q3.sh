#!/bin/bash
# B&B check of dnn.dscale on a solved sibling: tai25b (optimum 344,355,646 as the incumbent), tai50b's settings
# (experiments/estimates/tai50b/tai50b.par with dnn.depths = 1:23), without and with dnn.dscale = 0.5.
# The whole tree with a 14-minute limit each; if the baseline does not finish, the depth-2 subtree on the
# optimal path (-F 0:3,1:14) with a 7-minute limit each instead.
cd "$(dirname "$0")"
until [ "$(grep -c RESULT logs/tai50b_n7kid_w0.36_ds0.5.log 2>/dev/null)" = 2 ]; do sleep 5; done
OUT=$PWD/logs; PARD=$PWD
cd /home/user/wt-admm/latest
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
timeout 840 ./bnbdrv -p $PARD/bb_tai25b_base.par -u 344355646 ../qaplib/tai25b.dat > $OUT/bb_tai25b_base.log 2>&1
if [ $? = 124 ]; then
  timeout 420 ./bnbdrv -p $PARD/bb_tai25b_base.par -u 344355646 -F 0:3,1:14 ../qaplib/tai25b.dat > $OUT/bb_tai25b_F2_base.log 2>&1
  timeout 420 ./bnbdrv -p $PARD/bb_tai25b_ds0.5.par -u 344355646 -F 0:3,1:14 ../qaplib/tai25b.dat > $OUT/bb_tai25b_F2_ds0.5.log 2>&1
else
  timeout 840 ./bnbdrv -p $PARD/bb_tai25b_ds0.5.par -u 344355646 ../qaplib/tai25b.dat > $OUT/bb_tai25b_ds0.5.log 2>&1
fi
