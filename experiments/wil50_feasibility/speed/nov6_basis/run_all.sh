#!/bin/bash
# The measurement (README section 3): the warm set (run_warmset.sh, six families, 24
# children, production environment, AMX) in two passes of three configurations run
# at the same time, one per core, so that their CPU seconds are comparable:
#   pass 1: base (switch off, CPU 1), ship (DNN_SHIP_BASIS=1: one filter step after
#           the seed, CPU 2), ship0 (DNN_SHIP_BASIS_STEPS=0: the seed's Ritz pairs as
#           they are, CPU 3)
#   pass 2: base2 and ship2 (replicates, cores swapped: CPU 2, 1) and ship16 (the
#           parent's basis rounded to bf16 first, dnntest -Wb 16, as a state file with
#           QAP_DNN_SHIP_BASIS_BITS=16 would carry it; CPU 3)
# BIN: a copy of the dnntest built from the code commit (so that a rebuild during
# the run changes nothing). Logs: warm/<config>_<family>.log; then summarize.py.
D=$(cd "$(dirname "$0")" && pwd)
BIN=${BIN:?a copy of the dnntest to measure}
export BIN
FAMS=${FAMILIES:-"f9 f6 p10 p13 p2 p7"}
pass() {
  FAMILIES="$FAMS" CPU=$2 CONFIGS="$1" bash "$D/run_warmset.sh" &
  FAMILIES="$FAMS" CPU=$4 CONFIGS="$3" bash "$D/run_warmset.sh" &
  FAMILIES="$FAMS" CPU=$6 CONFIGS="$5" EXTRA="$7" bash "$D/run_warmset.sh" &
  wait
}
echo "$(date -u +%FT%TZ) pass 1"
pass "base@DNN_SHIP_BASIS=0" 1 "ship@DNN_SHIP_BASIS=1" 2 "ship0@DNN_SHIP_BASIS=1+DNN_SHIP_BASIS_STEPS=0" 3 ""
echo "$(date -u +%FT%TZ) pass 2"
pass "base2@DNN_SHIP_BASIS=0" 2 "ship2@DNN_SHIP_BASIS=1" 1 "ship16@DNN_SHIP_BASIS=1" 3 "-Wb 16"
echo "$(date -u +%FT%TZ) done"
# each configuration against the base of its own pass (same time, other cores)
python3 "$D/summarize.py" base ship ship0 --table > "$D/warm/summary_pass1.txt"
python3 "$D/summarize.py" base2 ship2 ship16 --table > "$D/warm/summary_pass2.txt"
# the replicates across passes (same configuration, other core, other time): the noise
{ python3 "$D/summarize.py" base base2; python3 "$D/summarize.py" ship ship2; } > "$D/warm/summary_noise.txt"
