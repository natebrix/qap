#!/bin/zsh
cd ~/dev/qap/latest
export VECLIB_MAXIMUM_THREADS=1 QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_BETA=0.21 QAP_DNN_REFRESH=120 QAP_DNN_WARM_Y=2 QAP_DNN_EXT=20 QAP_DNN_GU_ZONE=0.003 QAP_DNN_SHIP=32 DNN_ADMM_CROSS=2 DNN_CORE_LAZYCERT=1 DNN_BETA_SCHED=w=0.36
D=$HOME/dev/qap/experiments/tho40_estimate/depth2
c=$1
QAP_DNN_TRACE=$D/d2_26_$c.trace ./bnbdrv -p PARAM/nug_deep_xbar9_all.par -x 1:7:1:200:1 -u 240516 -F 26:$c -d 2 ../qaplib/tho40t.dat > $D/d2_26_$c.log 2>&1
