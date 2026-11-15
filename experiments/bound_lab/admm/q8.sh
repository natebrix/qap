#!/bin/bash
# queue 8 (tho40t, same parent and 36 children): w = 0.15 with dnn.cross_first = 5 under the production checks;
# w = 0.1 and 0.07 with a certificate at every iteration; and a give-up node warm for 200 iterations at
# w = 0.36 and 0.15 (the parent 26:10,8:0,39:32,38:7 itself as the warm child of 26:10,8:0,39:32, 600 cold)
cd "$(dirname "$0")"
while pgrep -f "q7.sh|bb3.sh" > /dev/null; do sleep 5; done
D=../../../latest/dnntest
P="-x 26:10,8:0,39:32,38:7 -J 600 -Q b -H 200 -i 200 -u 240516 -T 1 -v"
$D -p tho40_w0.15_cf5.par $P -k 10 -K "37:39;36:39;35:39;34:39;33:39;32:39;31:39;30:39;29:39;28:39;27:39;25:39;24:39;23:39;22:39;21:39;20:39;19:39;18:39;17:39;16:39;15:39;14:39;13:39;12:39;11:39;10:39;9:39;7:39;6:39;5:39;4:39;3:39;2:39;1:39;0:39" ../../../qaplib/tho40t.dat > logs/tho40_d4_w0.15_cf5.log 2>&1
$D -p tho40_w0.1.par $P -k 1 -K "37:39;36:39;35:39;34:39;33:39;32:39;31:39;30:39;29:39;28:39;27:39;25:39;24:39;23:39;22:39;21:39;20:39;19:39;18:39;17:39;16:39;15:39;14:39;13:39;12:39;11:39;10:39;9:39;7:39;6:39;5:39;4:39;3:39;2:39;1:39;0:39" ../../../qaplib/tho40t.dat > logs/tho40_d4_w0.1_k1.log 2>&1
for w in base w0.15; do
  $D -p tho40_$w.par -x 26:10,8:0,39:32 -J 600 -Q b -H 0 -i 200 -u 240516 -T 1 -v -k 10 -K "38:7" ../../../qaplib/tho40t.dat > logs/tho40_d3kid_$w.log 2>&1
done
$D -p tho40_w0.07.par $P -k 1 -K "37:39;36:39;35:39;34:39;33:39;32:39;31:39;30:39;29:39;28:39;27:39;25:39;24:39;23:39;22:39;21:39;20:39;19:39;18:39;17:39;16:39;15:39;14:39;13:39;12:39;11:39;10:39;9:39;7:39;6:39;5:39;4:39;3:39;2:39;1:39;0:39" ../../../qaplib/tho40t.dat > logs/tho40_d4_w0.07_k1.log 2>&1
