#!/bin/bash
# The default path is unchanged: the dnntest built before this study (OLD, from
# the branch point a4ad1db) and the new one, DNN_SHIP_BASIS unset, must print the
# same trace (every "it" line and RESULT, CPU fields stripped) on nug20's root,
# on wil50 all32 node 19 (cold) and on one warm child through -K ... -Q b (f9's
# parent and its child 37:21, the search's parent-Y start), each under two
# environments: the plain one of ../nov6_core/check_identity.sh and the production
# one (bestcore's DNN switches). Then, with the switch on, the parent's own trace
# of the warm-child run must still equal the old binary's (a cold call is never
# seeded; only its saved state grows by the basis).
# Run from anywhere: OLD=/path/to/old/dnntest check_identity.sh  (CPUA, CPUB: 1, 2;
# CPUA = CPUB runs the old and new binaries one after the other on that core)
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
OLD=${OLD:?path of the old dnntest}
NEW=${NEW:-$R/latest/dnntest}
CPUA=${CPUA:-1}; CPUB=${CPUB:-2}
BASE="OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1 DNN_AMX=1 DNN_R32_NONCERT=1"
PROD="$BASE DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 DNN_ADMM_CROSS=2 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1"
F="-T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40 -H 300 -v"
L=$D/logs/ident
mkdir -p "$L"
cd "$R" || exit 1
strip() { grep -E '^(it |RESULT)' "$1" | sed -E 's/ cpu [0-9.]+//; s/ cpu=[0-9.]+( ms_per_iter=[0-9.]+)?//; s/ seed=-?[0-9]+//'; }
# parent part of a -K run: everything before the first CHILD-SPEC
parent() { sed '/^CHILD-SPEC/,$d' "$1" | grep -E '^(it |RESULT)' | sed -E 's/ cpu [0-9.]+//; s/ cpu=[0-9.]+ ms_per_iter=[0-9.]+//'; }
res() { grep ^RESULT "$1" | sed -E 's/.*node=([0-9]+).*lb=([^ ]+).*status=([^ ]+) iters=([0-9]+) cpu=([0-9.]+).*/node \1 lb \2 \3 iters \4 cpu \5/'; }
# name|fixes|extra args (instance wil50 unless nug20)
for spec in "nug20|-|-i 200 -u 2570 qaplib/nug20.dat" \
            "w19|26:27,8:38,34:28,37:21|-i 2500 -u 48816 qaplib/wil50.dat" \
            "f9kid|26:27,8:38,34:28|-i 1500 -u 48816 -K 37:21 -Q b qaplib/wil50.dat"; do
  IFS='|' read -r name x args <<< "$spec"
  xarg=(); [ "$x" != "-" ] && xarg=(-x "$x")
  for envn in base prod; do
    E=$BASE; [ $envn = prod ] && E=$PROD
    if [ "$CPUA" = "$CPUB" ]; then
      env $E taskset -c $CPUA "$OLD" $F "${xarg[@]}" $args > "$L/${name}_${envn}_old.log" 2>&1
    else
      env $E taskset -c $CPUA "$OLD" $F "${xarg[@]}" $args > "$L/${name}_${envn}_old.log" 2>&1 &
    fi
    env $E taskset -c $CPUB "$NEW" $F "${xarg[@]}" $args > "$L/${name}_${envn}_new.log" 2>&1
    wait
    a=$L/${name}_${envn}_old.log b=$L/${name}_${envn}_new.log
    if cmp -s <(strip "$a") <(strip "$b"); then echo "$name/$envn: new identical to old ($(strip "$b" | wc -l) trace lines, CPU stripped)"
    else echo "$name/$envn: new DIFFERS from old"; diff <(strip "$a") <(strip "$b") | head -6; fi
    for v in old new; do res "$L/${name}_${envn}_$v.log" | sed "s/^/  $v: /"; done
  done
done
# the switch on: the cold parent of the warm-child run is unchanged
env $PROD DNN_SHIP_BASIS=1 taskset -c $CPUB "$NEW" $F -x 26:27,8:38,34:28 -i 1500 -u 48816 -K 37:21 -Q b qaplib/wil50.dat > "$L/f9kid_prod_on.log" 2>&1
if cmp -s <(parent "$L/f9kid_prod_old.log") <(parent "$L/f9kid_prod_on.log"); then echo "f9kid/prod, switch on: the parent's trace identical to old ($(parent "$L/f9kid_prod_on.log" | wc -l) lines)"
else echo "f9kid/prod, switch on: the parent's trace DIFFERS"; fi
res "$L/f9kid_prod_on.log" | sed "s/^/  on: /"
grep -h "^basis seed" "$L/f9kid_prod_on.log" | sed 's/^/  /'
