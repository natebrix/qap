#!/bin/sh
# small.sh INSTANCE GRID ITERS "options": the phase on a small instance's root (seconds), bound at marks
W=/home/user/wt-cpfast
inst=$1; grid=$2; iters=$3; shift 3
${BIN:-$W/bin/dnncuts.v2} ${PAR:+-p $PAR} -phase -x ${SPEC:--} -c ${CUTOFF:-100000000} -grid $grid -side 1 -I ${FIRST:-100} -k 500 -warmY -cpsep 50:500 -cpiters $iters $* -T 1 $W/qaplib/$inst.dat | \
 awk -v marks="${MARKS:-50 100 200 400 800 1600 3200 6400}" 'BEGIN{n=split(marks,M," ")} /^first solve/{printf "first %s ",$4} /^it /{for(i=1;i<=n;i++) if($2==M[i]) printf "%s:%s ",$2,$4} /^extrap guard/{g++; if($NF=="REVERT") r++} /^FINAL/{printf "| jumps %d reverts %d cpu %s\n", g, r, $NF}'
