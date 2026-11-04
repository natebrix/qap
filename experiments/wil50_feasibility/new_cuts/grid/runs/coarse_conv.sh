#!/bin/bash
# converged (1,500 it.) plain DNN of node 13's mid5 coarse problem, its lowest
# two coarse children for facility 40 (cells 24 and 34), sequentially on core 2
cd /home/user/qap
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_AMX=1 DNN_CHOLQR_PASSES=1
G=experiments/wil50_feasibility/new_cuts/grid
for spec in 21:28,8:41,28:19 21:28,8:41,28:19,24:40 21:28,8:41,28:19,34:40; do
  t=${spec//[:,]/_}; f=$G/runs/coarseconv_$t.log
  grep -q RESULT $f 2>/dev/null && continue
  taskset -c 2 latest/dnntest -x $spec -c 1e12 -T 1 -d 3 -L 2 -r 40 -i 1500 -H 0 $G/inst/wil50_mid5only.dat > $f.tmp 2>&1 && mv $f.tmp $f
done
