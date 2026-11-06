#!/bin/bash
# The fallback: a shipped block that does not fit is rejected and the child refreshes
# as without the switch. nug30's root (300 iterations, no give-up) as the parent, two
# children through -K ... -Q b, production environment, DNN_SHIP_BASIS=1, with the
# block intact, replaced by pseudo-random vectors (DNN_SHIP_BASIS_FAULT=1) and with
# its coordinates reversed (=2, a block in the wrong basis); then the switch off.
# The faulty runs must report "rejected" and give the switch-off children's trace.
# Logs: logs/fault/. CPU: the core (1).
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
BIN=${BIN:-$R/latest/dnntest}
CPU=${CPU:-1}
L=$D/logs/fault
mkdir -p "$L"
PROD="OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1 DNN_AMX=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 DNN_ADMM_CROSS=2 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1 DNN_R32_NONCERT=1"
F="-T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40 -H 300 -i 1500 -J 300 -v -v -u 6124 -K 3:5,7:7 -Q b $R/qaplib/nug30.dat"
for v in "ok@DNN_SHIP_BASIS=1" "random@DNN_SHIP_BASIS=1 DNN_SHIP_BASIS_FAULT=1" "reversed@DNN_SHIP_BASIS=1 DNN_SHIP_BASIS_FAULT=2" "off@DNN_SHIP_BASIS=0"; do
  name=${v%%@*}; envs=${v#*@}
  env $PROD $envs taskset -c $CPU "$BIN" $F > "$L/$name.log" 2>&1
  echo "$name: $(grep -h '^basis seed' "$L/$name.log" | sed -E 's/basis seed (\S+): kb ([0-9]+) nneg ([0-9]+) res (\S+) .*first R-step (\S+) ms/\1 (kb \2, res \4, first R-step \5 ms)/' | tr '\n' ';')"
  grep ^RESULT "$L/$name.log" | sed -E 's/.*node=([0-9]+).*lb=([^ ]+).*status=([^ ]+) iters=([0-9]+) cpu=([0-9.]+).*/  node \1 lb \2 \3 iters \4 cpu \5/'
done
kids() { sed -n '/^CHILD-SPEC/,$p' "$1" | grep -E '^(it |RESULT)' | sed -E 's/ cpu [0-9.]+//; s/ cpu=[0-9.]+//; s/ seed=-?[0-9]+//'; }
for name in random reversed; do
  if cmp -s <(kids "$L/$name.log") <(kids "$L/off.log"); then echo "$name: the children's trace equals the switch-off run"; else echo "$name: the children's trace DIFFERS from the switch-off run"; fi
done
