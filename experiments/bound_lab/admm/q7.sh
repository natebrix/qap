#!/bin/bash
# queue 7: tho40t parent and 36 children as q2.sh, a certificate at every iteration, the warm-call beta factor
# w = 0.15, 0.25, 0.5 (production 0.36); then the tai30b subtree comparison (bb3.sh's body)
cd "$(dirname "$0")"
while pgrep -f "q2.sh|q6.sh" > /dev/null; do sleep 5; done
D=../../../latest/dnntest
P="-x 26:10,8:0,39:32,38:7 -J 600 -Q b -H 200 -i 200 -u 240516 -T 1 -v"
for w in 0.25 0.15 0.5; do
  $D -p tho40_w$w.par $P -k 1 -K "37:39;36:39;35:39;34:39;33:39;32:39;31:39;30:39;29:39;28:39;27:39;25:39;24:39;23:39;22:39;21:39;20:39;19:39;18:39;17:39;16:39;15:39;14:39;13:39;12:39;11:39;10:39;9:39;7:39;6:39;5:39;4:39;3:39;2:39;1:39;0:39" ../../../qaplib/tho40t.dat > logs/tho40_d4_w${w}_k1.log 2>&1
done
sed -i 's/^sleep 30$//; s/^while pgrep.*$//' bb3.sh
./bb3.sh
