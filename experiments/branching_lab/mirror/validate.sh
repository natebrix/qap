#!/bin/bash
# the small searches validate.py checks (bnbdrv.hook, QAP_MIRROR_CHECK=1), one at a time
H=$(cd "$(dirname "$0")" && pwd); Q=/home/user/wt-cognac/qaplib; export QAP_MIRROR_CHECK=1
run() { NAME=$1; INST=$2; UB=$3; export QAP_MIRROR="$4"; STORE=$H/runs/val_store $H/run.sh $NAME hook -p $H/small.par -u $UB $Q/$INST.dat; }
for mu in 0.25 1 8; do run val_tai10b_mu$mu tai10b 1183761 "0,1,$mu;2,3,$mu;4,5,$mu;0,5,$mu;1,7,$mu"; done
for mu in 0.5 4; do run val_nug08_mu$mu nug08 216 "0,1,$mu;2,5,$mu;3,6,$mu;1,4,$mu"; done
for mu in 0.5 2; do run val_tai12b_mu$mu tai12b 39464926 "0,1,$mu;2,3,$mu;4,7,$mu;5,9,$mu;1,10,$mu"; done
