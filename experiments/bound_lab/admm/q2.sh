#!/bin/bash
# part B: tho40t give-up parent 26:10,8:0,39:32,38:7 (depth 4, m = 36; production: 406 short, 36 children all fathom)
# parent cold for 600 iterations (-J 600), then its 36 children warm (Z and conditioned Y: -Q b), give-up 200:
#   k1: a certificate at every iteration (the first iteration at which each child's certificate crosses)
#   prod: tho40.par as is (checks every 10, dnn.cross = 2); cf5 / cf3: with dnn.cross_first = 5 / 3
cd "$(dirname "$0")"
while pgrep -f "q1.sh|q3.sh" > /dev/null; do sleep 5; done
D=../../../latest/dnntest
P="-x 26:10,8:0,39:32,38:7 -J 600 -Q b -H 200 -i 200 -u 240516 -T 1 -v"
$D -p tho40_base.par $P -k 1 -K "37:39;36:39;35:39;34:39;33:39;32:39;31:39;30:39;29:39;28:39;27:39;25:39;24:39;23:39;22:39;21:39;20:39;19:39;18:39;17:39;16:39;15:39;14:39;13:39;12:39;11:39;10:39;9:39;7:39;6:39;5:39;4:39;3:39;2:39;1:39;0:39" ../../../qaplib/tho40t.dat > logs/tho40_d4_k1.log 2>&1
$D -p tho40_base.par $P -k 10 -K "37:39;36:39;35:39;34:39;33:39;32:39;31:39;30:39;29:39;28:39;27:39;25:39;24:39;23:39;22:39;21:39;20:39;19:39;18:39;17:39;16:39;15:39;14:39;13:39;12:39;11:39;10:39;9:39;7:39;6:39;5:39;4:39;3:39;2:39;1:39;0:39" ../../../qaplib/tho40t.dat > logs/tho40_d4_prod.log 2>&1
$D -p tho40_cf5.par $P -k 10 -K "37:39;36:39;35:39;34:39;33:39;32:39;31:39;30:39;29:39;28:39;27:39;25:39;24:39;23:39;22:39;21:39;20:39;19:39;18:39;17:39;16:39;15:39;14:39;13:39;12:39;11:39;10:39;9:39;7:39;6:39;5:39;4:39;3:39;2:39;1:39;0:39" ../../../qaplib/tho40t.dat > logs/tho40_d4_cf5.log 2>&1
$D -p tho40_cf3.par $P -k 10 -K "37:39;36:39;35:39;34:39;33:39;32:39;31:39;30:39;29:39;28:39;27:39;25:39;24:39;23:39;22:39;21:39;20:39;19:39;18:39;17:39;16:39;15:39;14:39;13:39;12:39;11:39;10:39;9:39;7:39;6:39;5:39;4:39;3:39;2:39;1:39;0:39" ../../../qaplib/tho40t.dat > logs/tho40_d4_cf3.log 2>&1
