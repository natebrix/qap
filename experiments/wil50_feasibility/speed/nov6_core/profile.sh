#!/bin/bash
# Profile of a production iteration at m = 47 (node 10, 39:27,30:48,25:49) on CPU 3:
# cold without AMX, cold with AMX, cold with the production `fast` switches,
# warm-started child (parent 39:27,30:48 at m = 48, child 25:49 from Z + conditioned Y),
# and a depth-4 fathoming node (all32 line 19). dnntest -v -v prints the phase
# profile, the subspace step split and the form-R split. One process at a time.
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
BIN=${BIN:-$R/latest/dnntest}
CPU=${CPU:-3}
FLAGS="-T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40 -H 300 -u 48816 -v -v"
BASE="OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1"
PROD="DNN_AMX=1 DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.7"
run() { # name env... -- args
  local name=$1; shift
  local envs=()
  while [ "$1" != "--" ]; do envs+=("$1"); shift; done; shift
  echo "$(date -u +%T) $name: load $(cut -d' ' -f1-3 /proc/loadavg)" >> "$D/logs/profile_index.txt"
  env $BASE "${envs[@]}" taskset -c $CPU "$BIN" $FLAGS "$@" "$R/qaplib/wil50.dat" > "$D/logs/prof_$name.log" 2>&1
}
cd "$R" || exit 1
run cold_sgemm        -- -i 300 -x 39:27,30:48,25:49
run cold_amx  DNN_AMX=1 -- -i 300 -x 39:27,30:48,25:49
run cold_prod $PROD   -- -i 300 -x 39:27,30:48,25:49
run warm_prod $PROD   -- -i 1500 -J 300 -K 25:49 -Q b -x 39:27,30:48
run d4n19_prod $PROD  -- -i 1500 -x 26:27,8:38,34:28,37:21
echo "$(date -u +%T) done" >> "$D/logs/profile_index.txt"
