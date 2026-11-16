#!/bin/bash
# second resume of the tai30b search (bnb_tai30b_cut.state2 -> state3), 420 CPU-s
L=/home/user/qap/experiments/bound_lab/regions/logs
while pgrep -f "bnbdrv -p /home/user/qap/experiments/bound_lab/regions" >/dev/null; do sleep 3; done
[ -f $L/bnb_tai30b_cut.state2 ] || exit 0
cd /home/user/wt-regions/latest
QAP_NODE_TRACE=$L/bnb_tai30b_cut.trace3 ./bnbdrv -p /home/user/qap/experiments/bound_lab/regions/bnb/tai30b_cut.par -r $L/bnb_tai30b_cut.state2 -o $L/bnb_tai30b_cut.state3 -t 420 ../qaplib/tai30b.dat > $L/bnb_tai30b_cut.log3 2>&1
