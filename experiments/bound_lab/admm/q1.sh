#!/bin/bash
# queue 1: tai30b root sweep (p, beta), then tai50b depth-3 nodes 7 and 3 (NOTES.md section 2), cold, 600 iterations
cd "$(dirname "$0")"
WT=/home/user/wt-admm/latest/dnntest
Q=../../../qaplib
run() { # tag par inst ub iters extra...
  tag=$1; par=$2; inst=$3; ub=$4; it=$5; shift 5
  [ -s logs/$tag.log ] && grep -q RESULT logs/$tag.log && return
  $WT -p $par -T 1 -H 0 -i $it -k 10 -u $ub -v "$@" $Q/$inst.dat > logs/$tag.log 2>&1
}
run tai30b_root_ds0.5 par_ds0.5.par tai30b 637117113 1500
run tai30b_root_ds0.75 par_ds0.75.par tai30b 637117113 1500
run tai30b_root_ds0.25 par_ds0.25.par tai30b 637117113 1500
run tai30b_root_ds1.0 par_ds1.0.par tai30b 637117113 1500
run tai30b_root_ds0.5_b0.35 par_ds0.5.par tai30b 637117113 1500 -B 0.35
run tai30b_root_ds0.5_b0.12 par_ds0.5.par tai30b 637117113 1500 -B 0.12
run tai50b_n7_base base.par tai50b 458821517 600 -x 0:44,49:49,25:31
run tai50b_n7_ds0.5 par_ds0.5.par tai50b 458821517 600 -x 0:44,49:49,25:31
run tai50b_n3_base base.par tai50b 458821517 600 -x 27:0,44:23,6:19
run tai50b_n3_ds0.5 par_ds0.5.par tai50b 458821517 600 -x 27:0,44:23,6:19
