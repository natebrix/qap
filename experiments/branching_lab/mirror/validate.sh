#!/bin/bash
# the small searches validate.py checks (bnbdrv.hook, QAP_MIRROR_CHECK=1), one at a time; each from optimum + 1
# (+ 2 on nug08, whose values are even) and from a loose incumbent (about 3% above), so more nodes are bounded
H=$(cd "$(dirname "$0")" && pwd); Q=/home/user/wt-cognac/qaplib; export QAP_MIRROR_CHECK=1
run() { NAME=$1; INST=$2; UB=$3; export QAP_MIRROR="$4"; STORE=$H/runs/store_val $H/run.sh $NAME hook -p $H/small.par -u $UB $Q/$INST.dat; }
for ub in 1183761 1220000; do for mu in 0.25 1 8; do run val_tai10b_u${ub}_mu$mu tai10b $ub "0,1,$mu;2,3,$mu;4,5,$mu;0,5,$mu;1,7,$mu;6,8,$mu"; done; done
for ub in 216 222; do for mu in 0.5 4; do run val_nug08_u${ub}_mu$mu nug08 $ub "0,1,$mu;2,5,$mu;3,6,$mu;1,4,$mu;0,7,$mu"; done; done
for ub in 39464926 40700000; do for mu in 0.5 2; do run val_tai12b_u${ub}_mu$mu tai12b $ub "0,1,$mu;2,3,$mu;4,7,$mu;5,9,$mu;1,10,$mu;6,11,$mu"; done; done
