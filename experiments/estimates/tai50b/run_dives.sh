#!/bin/bash
# tai50b: Knuth dives (bnbdrv -e) with the current settings (tai50b.par, settings.sh),
# one process per core with distinct seeds, for HOURS of wall, then the pooled estimate.
#
#   experiments/estimates/tai50b/run_dives.sh [HOURS] [JOBS] [FIRST_SEED]
#   experiments/estimates/tai50b/run_dives.sh 8          # overnight on the performance cores
#
# estimate.par (the estimator reads it from its working directory): dives start at level 1
# (the root's children chosen by score), 10,000 dives at most, the per-level table printed
# after every dive, score exponent 1.5, as in the Nov 8 tuned run
# (../../tai50b_feasibility/NOTES.md section 5). The dive sampler is in the settings
# (est.xbar = 1, est.xbar_mix = 0.3). A process stopped at the wall limit keeps every
# finished dive. Seeds FIRST_SEED.. (default 101..), so a second night adds new dives:
# run it again with another FIRST_SEED and pool all directories.
# Output: experiments/estimates/tai50b/dives/sSEED/{est.out,dives.txt}; pool.txt and the
# CSVs of ../../tai50b_feasibility/estimate/pool.py in dives/.
set -u
HOURS=${1:-8}
R=$(cd "$(dirname "$0")/../../.." && pwd); cd "$R" || exit 1
. experiments/estimates/tai50b/settings.sh || exit 1
fullrun_check || exit 1
if [ "$(uname -s)" = Darwin ]; then
  P=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.physicalcpu)
  command -v caffeinate > /dev/null && { caffeinate -ims -w $$ & }
else P=$(nproc); fi
JOBS=${2:-$P}; S0=${3:-101}
D=$R/experiments/estimates/tai50b/dives; mkdir -p "$D"
SECS=$(python3 -c "print(int($HOURS*3600))")
seeds=()
for k in $(seq 0 $((JOBS - 1))); do
  s=$((S0 + k)); seeds+=("$s"); d=$D/s$s
  [ -f "$d/est.out" ] && { echo "seed $s: $d/est.out exists; skipped (pick another FIRST_SEED)"; continue; }
  mkdir -p "$d" && cp experiments/estimates/tai50b/estimate.par "$d/"
  # the machine, since CPU seconds differ by machine (node counts pool across machines, CPU only within one)
  { if [ "$(uname -s)" = Darwin ]; then sysctl -n machdep.cpu.brand_string; else grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | sed 's/^ //'; fi
    echo "$(date -u +%FT%TZ) $(git rev-parse --short HEAD) ${HOURS} h"; } > "$d/host.txt"
  ( cd "$d" && QAP_EST_SEED=$s QAP_EST_HIST=$d/dives.txt timeout -s TERM "$SECS" "$R/latest/bnbdrv" -e -p "$FULLRUN_PAR" \
      -u "$FULLRUN_UB" "$FULLRUN_INSTANCE" > "$d/est.out" 2>&1 ) &
done
echo "$(date '+%F %T') tai50b: ${#seeds[@]} estimator processes (seeds ${seeds[*]}), ${HOURS} h; output in $D"
wait
all=$(cd "$D" && ls -d s* | sed 's/^s//' | sort -n | tr '\n' ' ')
python3 -I experiments/tai50b_feasibility/estimate/pool.py "$D" $all | tee "$D/pool.txt"
echo; echo "To share (est.out and dives.txt are ignored while running): git add -f experiments/estimates/tai50b/dives && git commit -m 'tai50b dives' && git push"
