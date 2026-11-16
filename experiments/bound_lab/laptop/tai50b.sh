#!/bin/bash
# tai50b whole, through the coordinator's local mode, with the cluster-region cut phase
# (bound_lab/README.md section 1). A rehearsal, not a run of record: a build with today's
# patches (cut.regions, the pooled-solve retry of tai_solves/pool_retry.patch, cut.giveup_depth),
# settings nobody has run on tai50b in this exact form (below), one machine.
#
# From the repository's root (overnight.sh starts it):
#   caffeinate -i experiments/bound_lab/laptop/tai50b.sh
# Options (environment): INST (tai50b; tai35b and tai40b are the tests), UB (best known + 1 from
# QAPSOLN/, so the search must find the best known value itself), SLOTS (all cores), OUT
# (runs_INST/ here), MAX_WALL_H (20), MIN_FREE_GB (10: the run is stopped when the disk has less).
# A run whose results are there is not repeated. An interrupted run starts over (the coordinator's
# store is kept in OUT/store for a look, and removed by the next start).
#
# Settings: bound_lab/tai50b_size/par/A_wc.par (tai50b's settings, the region file, a phase at every
# DNN miss, pools inherited, dnn.beta_sched = w=0.36,wc=2.88) with four changes, so that one file
# serves the root and the rest:
#   cut.stall = 0, cut.cp_iters = 6000   the root's phase runs long (6,000 iterations), as the
#                                        roots of tai35b, tai40b and tai50b did
#   cut.giveup = 2400:0, cut.giveup_depth = 1   below the root a phase stops at 2,400 iterations;
#                                        the cloud runs stopped stalled phases earlier (cut.stall = 3)
#   ship.bits = 32                       parents' DNN states in fp32, half the disk
set -uo pipefail
H=$(cd "$(dirname "$0")" && pwd)
cd "$H/../../.."
R=$(pwd)
INST=${INST:-tai50b}
SLOTS=${SLOTS:-$( (sysctl -n hw.ncpu 2>/dev/null || nproc) )}
OUT=${OUT:-$H/runs_$INST}
MAX_WALL_H=${MAX_WALL_H:-20}
MIN_FREE_GB=${MIN_FREE_GB:-10}
PY=${PYTHON:-python3}
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
REG=$R/experiments/bound_lab/regions/regions/$INST.comb.txt
[ -f "$REG" ] || REG=$R/experiments/bound_lab/tai_solves/regions/$INST.comb.txt
[ -f "$REG" ] || { echo "no region file for $INST"; exit 1; }
BEST=$(awk 'NR==1{gsub(/[^0-9]/, "", $2); print $2}' "$R/QAPSOLN/$INST.sln")
UB=${UB:-$(( BEST + 1 ))}
N=$(awk 'NF{gsub(/[^0-9]/, "", $1); print $1; exit}' "$R/qaplib/$INST.dat")
RUN=$INST-lab
mkdir -p "$OUT"
say() { echo "$(date '+%m-%d %H:%M:%S') $*" | tee -a "$OUT/run.log"; }
free_gb() { df -Pk "$OUT" | awk 'NR==2{printf "%d", $4 / 1048576}'; }

res=$OUT/results.txt
if [ -s "$res" ]; then say "$INST: done before ($res)"; exit 0; fi
say "$INST: n = $N, best known $BEST, incumbent $UB, $SLOTS slots, $(free_gb) GB free, commit $(git rev-parse --short HEAD)"
if [ "$(free_gb)" -lt $(( MIN_FREE_GB * 2 )) ]; then say "less than $(( MIN_FREE_GB * 2 )) GB free: not started"; exit 1; fi
make solver PYTHON="$PY" > "$OUT/build.log" 2>&1 || { say "build failed: $OUT/build.log"; exit 1; }

"$PY" - "$R/experiments/bound_lab/tai50b_size/par/A_wc.par" "$OUT/$INST.par" "$REG" "$N" << 'EOF'
import re, sys
base, out, reg, n = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
kv = {"dnn.depths": f"0:{n - 2}", "cut.regions": reg, "cut.stall": "0", "cut.cp_iters": "6000",
      "cut.giveup": "2400:0", "cut.giveup_depth": "1", "ship.bits": "32"}
t = open(base).read()
for k, v in kv.items():
    pat = rf"^{re.escape(k)} = .*$"
    t = re.sub(pat, lambda m: f"{k} = {v}", t, flags=re.M) if re.search(pat, t, re.M) else t.rstrip("\n") + f"\n{k} = {v}\n"
open(out, "w").write(t)
EOF

d=$OUT/store; rm -rf "$d"
# the disk watch: stop the run when the disk runs low (one saved state per branched node)
( while sleep 120; do
    if [ "$(free_gb)" -lt "$MIN_FREE_GB" ]; then
      say "disk below $MIN_FREE_GB GB: stopping the run"
      "$PY" -m coord.qapctl stop "$RUN" --backend "local:$d" >> "$OUT/coord.log" 2>&1
      break
    fi
  done ) &
watch=$!
trap 'kill $watch 2> /dev/null' EXIT

say "$INST: running (the root alone takes 10 to 70 minutes on one core; status: $OUT/coord.log)"
if "$PY" -m coord.qapctl local "$RUN" --backend "local:$d" --instance "$R/qaplib/$INST.dat" \
     --params "$OUT/$INST.par" --ub "$UB" --slots "$SLOTS" \
     --budget 600 --target-unit-cpu 600 --steal-after 60 --checkpoint-every 120 --ckpt-upload 300 \
     --max-wall $(( MAX_WALL_H * 3600 )) \
     --trace --visibility private --workdir "$d/work" > "$OUT/coord.log" 2>&1 \
   && [ -s "$d/store/runs/$RUN/results.txt" ]; then
  cp "$d/store/runs/$RUN/results.txt" "$res"
  "$PY" -m coord.qapctl export "$RUN" --backend "local:$d" --out "$OUT/export" >> "$OUT/coord.log" 2>&1 \
    && rm -rf "$d" || say "export failed, the store is kept in $d"
  say "$INST: $(grep -E 'Nodes:|coordinator:' "$res" | tr -s ' ' | tr '\n' ' ' | cut -c1-240)"
  best=$(grep -oE 'best objective [0-9]+' "$res" | tail -1 | awk '{print $3}')
  if [ "$best" = "$BEST" ]; then
    say "$INST: ended at the best known value $BEST with the tree closed (a rehearsal: see the header)"
  else
    say "$INST: ended at $best, NOT the best known value $BEST: look at $res before anything else"
  fi
else
  say "$INST: did not finish (stopped, out of time or failed): see $OUT/coord.log; the store is in $d"
  exit 1
fi
