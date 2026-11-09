#!/bin/bash
# Knuth estimate of the tai50b tree (bnbdrv -e) with the nug30 record recipe
# (experiments/nug30_solve/README.md, Nov 7: rule 9 through depth 6, DNN at
# depths 1-7, give-up 200) and the DNN-marginal dive sampler with a 30% uniform
# share. Dives from level 1; the per-level table is printed after every dive
# (estimate.par: 1 10000 1 1 1.5 1000), so per-level variances can be computed.
# Usage: [X=D1:D2:G:H:W TAG=name] run_est.sh SEED   (one process per core, distinct seeds)
# X: the -x argument (default 1:7:1:200:1, the nug30 record recipe); TAG names the run directory.
seed=$1; X=${X:-1:7:1:200:1}; TAG=${TAG:-}; here=$(cd "$(dirname "$0")" && pwd); L=$here/../../../latest
d=$here/${TAG:+$TAG/}s$seed; mkdir -p $d && cd $d && ln -sfn $L/PARAM PARAM && ln -sfn $L/SYM SYM && cp $here/estimate.par .
export OPENBLAS_NUM_THREADS=1
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 \
       QAP_DNN_BETA=0.21 QAP_DNN_REFRESH=120 QAP_DNN_WARM_Y=2 DNN_CHOLQR_PASSES=1 \
       QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003 QAP_DNN_SHIP=32 DNN_R32_NONCERT=1 \
       DNN_ADMM_CROSS=2 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 \
       DNN_CORE_REFRESH_PARTIAL=1 DNN_BETA_SCHED=w=0.36
export QAP_EST_SEED=$seed QAP_EST_XBAR=1 QAP_EST_XBAR_MIX=0.3 QAP_EST_HIST=$d/dives.txt
exec $L/bnbdrv -e -p PARAM/nug_deep_xbar_all.par -x $X -u 458821517 $here/../../../qaplib/tai50b.dat > $d/est.out 2>&1
