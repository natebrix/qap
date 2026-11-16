mkpar () 
{ 
    cp $R/bnb/head.par par/$1.par;
    cat >> par/$1.par <<EOF
# bound_lab/regions/size50: experiments/estimates/tai50b/tai50b.par's values, then the cut block
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
dnn.beta_sched = w=0.36
dnn.depths = $4
dnn.maxgap = 1
dnn.giveup = $5
dnn.giveup_warm = $6
dnn.warm = 1
est.xbar = 1
est.xbar_mix = 0.3
# the cut phase on the cluster regions at every DNN miss, CP form, pools inherited
pool.on = 1
cut.regions = $R/regions/$3
cut.side = 0
cut.zone = 1
cut.warmy = 1
cut.rounds = 12
cut.iters = 200
cut.freeze = 0
cut.cp = 100:1000
cut.cp_iters = $7
cut.stall = $8
EOF

    shift 8;
    for l in "$@";
    do
        echo "$l" >> par/$1.par;
    done
}
