#!/bin/bash
# mkpar.sh NAME REGIONS DEPTHS CP_ITERS STALL [extra "key = value" lines...]: par/NAME.par = the rule block of
# ../regions/bnb/head.par + the tai50b sizing settings (../regions/size50/mkpar.sh, kidW288) with the phase
# penalty wc=2.88 from the start. REGIONS is a file in regions/. Later lines override earlier ones? No:
# a key may appear once, so extras must be new keys (use the sed in the callers to change a value).
D=$(cd $(dirname $0); pwd); N=$1; R=$2; DEPTHS=$3; CPI=$4; STALL=$5; shift 5
cp $D/../regions/bnb/head.par $D/par/$N.par
cat >> $D/par/$N.par <<EOT
# bound_lab/tai60b: experiments/estimates/tai50b/tai50b.par's values, then the cut block of regions/size50
dnn.gamma1 = 1.0
dnn.gamma = 1.0
dnn.cheb = 3
dnn.tight_lam = 2
dnn.beta = 0.21
dnn.refresh = 120
dnn.warm_y = 2
dnn.ext = 20
dnn.gu_zone = 0
ship.bits = 64
dnn.cross = 2
dnn.lazycert = 1
dnn.beta_sched = ${BS:-w=0.36,wc=2.88}
dnn.depths = $DEPTHS
dnn.maxgap = 1
dnn.giveup = 100
dnn.giveup_warm = 100
dnn.warm = 1
est.xbar = 1
est.xbar_mix = 0.3
# the cut phase on the cluster regions at every DNN miss, CP form, pools inherited
pool.on = 1
cut.regions = $D/regions/$R
cut.side = 0
cut.zone = ${ZONE:-1}
cut.warmy = 1
cut.rounds = 12
cut.iters = 200
cut.freeze = 0
cut.cp = 100:1000
cut.cp_iters = $CPI
cut.stall = $STALL
EOT
for l in "$@"; do echo "$l" >> $D/par/$N.par; done
