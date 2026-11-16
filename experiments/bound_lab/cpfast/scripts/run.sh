#!/bin/sh
# run.sh NAME SPEC [dnncuts options...]: one phase run of the worktree's dnncuts, log in logs/NAME.log
# (first line: the command). The first solve's state is cached in state/SPEC_I.warm.
D=/home/user/qap/experiments/bound_lab/cpfast
W=/home/user/wt-cpfast
name=$1; spec=$2; shift 2
tag=$(echo "$spec" | tr ':,' '__')
I=${FIRST:-300}
warm=$D/state/${tag}_I$I${FB:+_fb$FB}.warm
if [ -f "$warm" ]; then wopt="-loadwarm $warm"; else wopt="-savewarm $warm"; fi
cmd="${BIN:-$W/bin/dnncuts.v4} -p ${PAR:-$W/experiments/estimates/tho40/tho40.par} -phase -x $spec -c ${CUTOFF:-241500} -grid 5:8 -side 1 -I $I -k 1000 -warmY ${FB:+-firstbeta $FB} $wopt -savepool $D/state/$name.pool $* -T 1 $W/qaplib/tho40t.dat"
echo "# $cmd" > $D/logs/$name.log
cd $W && $cmd >> $D/logs/$name.log 2>&1
echo "# exit $?" >> $D/logs/$name.log
