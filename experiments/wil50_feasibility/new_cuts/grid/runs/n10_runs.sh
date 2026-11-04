#!/bin/bash
# node 10: tri5 with dump, then objective splitting (mid5, alpha 0.5) on top of the triangles
cd /home/user/qap
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_AMX=1 DNN_CHOLQR_PASSES=1
G=experiments/wil50_feasibility/new_cuts/grid
grep -q FINAL $G/runs/n10_tri5.log 2>/dev/null || taskset -c 2 experiments/dnn_cuts/dnncuts -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 -i 300 -r 5 -delta 50 -nocap -freeze 0:600 -x 39:27,30:48,25:49 -dumpY $G/runs/n10_tri5 qaplib/wil50.dat > $G/runs/n10_tri5.log 2>&1
grep -q FINAL $G/runs/n10_tri5_mid5a1of2.log 2>/dev/null || taskset -c 2 experiments/dnn_cuts/dnncuts -grid 5:10 -c 1e12 -T 1 -k 1000 -v -I 400 -i 300 -r 5 -delta 100 -nocap -freeze 0:600 -x 39:27,30:48,25:49 $G/inst/wil50_mid5_a1of2.dat > $G/runs/n10_tri5_mid5a1of2.log 2>&1
