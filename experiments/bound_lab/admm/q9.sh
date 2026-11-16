#!/bin/bash
# queue 9 (tho40t): warm calls at 0.1 x base for the first 30 iterations, then 0.36 (DNN_ADMM_WARM_BETA = 0.1:30:0.36:1),
# with dnn.cross_first = 5: the 36 children under the production checks, and the give-up node warm for 200 iterations
cd "$(dirname "$0")"
while pgrep -f "bb4.sh" > /dev/null; do sleep 5; done
D=../../../latest/dnntest
DNN_ADMM_WARM_BETA=0.1:30:0.36:1 $D -p tho40_cf5.par -x 26:10,8:0,39:32,38:7 -J 600 -Q b -H 200 -i 200 -u 240516 -T 1 -v -k 10 -K "37:39;36:39;35:39;34:39;33:39;32:39;31:39;30:39;29:39;28:39;27:39;25:39;24:39;23:39;22:39;21:39;20:39;19:39;18:39;17:39;16:39;15:39;14:39;13:39;12:39;11:39;10:39;9:39;7:39;6:39;5:39;4:39;3:39;2:39;1:39;0:39" ../../../qaplib/tho40t.dat > logs/tho40_d4_wb0.1_30_cf5.log 2>&1
DNN_ADMM_WARM_BETA=0.1:30:0.36:1 $D -p tho40_cf5.par -x 26:10,8:0,39:32 -J 600 -Q b -H 0 -i 200 -u 240516 -T 1 -v -k 10 -K "38:7" ../../../qaplib/tho40t.dat > logs/tho40_d3kid_wb0.1_30.log 2>&1
