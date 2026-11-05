#!/bin/bash
# Validity checks with DNN_RPRSM set, on core 2, one job at a time.
#   usage: run_checks.sh RPRSM [bruteforce] [dnncheck]
# Logs: logs/brute_force_rprsm<RPRSM>.log, logs/dnncheck_rprsm<RPRSM>.log
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
RP=$1; shift
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
if [ "$RP" != 0 ]; then export DNN_RPRSM=$RP; else unset DNN_RPRSM; fi
cd "$R"
for what in "$@"; do
  case $what in
    bruteforce)
      taskset -c 2 python3 tests/brute_force.py --jobs 1 --solver-args "-x 0:20" \
        > "$D/logs/brute_force_rprsm$RP.log" 2>&1
      echo "bruteforce exit $?"; tail -n 3 "$D/logs/brute_force_rprsm$RP.log" ;;
    dnncheck)
      taskset -c 2 make dnncheck > "$D/logs/dnncheck_rprsm$RP.log" 2>&1
      echo "dnncheck exit $?"; tail -n 5 "$D/logs/dnncheck_rprsm$RP.log" ;;
  esac
done
