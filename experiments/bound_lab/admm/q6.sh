#!/bin/bash
# queue 6 (after the tho40 runs of q2.sh): the same tho40t parent and 36 children with dnn.dscale = 0.5
# (patched binary), every-iteration certificates and the production checks; then the patched binary with
# dnn.dscale = 0 as its own control for the production checks.
cd "$(dirname "$0")"
sleep 20
while pgrep -f "q2.sh|bb_q3.sh" > /dev/null; do sleep 5; done
WT=/home/user/wt-admm/latest/dnntest
P="-x 26:10,8:0,39:32,38:7 -J 600 -Q b -H 200 -i 200 -u 240516 -T 1 -v"
$WT -p tho40_ds0.5.par $P -k 1 -K "37:39;36:39;35:39;34:39;33:39;32:39;31:39;30:39;29:39;28:39;27:39;25:39;24:39;23:39;22:39;21:39;20:39;19:39;18:39;17:39;16:39;15:39;14:39;13:39;12:39;11:39;10:39;9:39;7:39;6:39;5:39;4:39;3:39;2:39;1:39;0:39" ../../../qaplib/tho40t.dat > logs/tho40_d4_ds0.5_k1.log 2>&1
$WT -p tho40_ds0.5.par $P -k 10 -K "37:39;36:39;35:39;34:39;33:39;32:39;31:39;30:39;29:39;28:39;27:39;25:39;24:39;23:39;22:39;21:39;20:39;19:39;18:39;17:39;16:39;15:39;14:39;13:39;12:39;11:39;10:39;9:39;7:39;6:39;5:39;4:39;3:39;2:39;1:39;0:39" ../../../qaplib/tho40t.dat > logs/tho40_d4_ds0.5_prod.log 2>&1
