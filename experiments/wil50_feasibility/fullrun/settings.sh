#!/bin/bash
# The full wil50 run's settings: one place for the launcher (not built yet),
# the sampling plan (roadmap/tree_size.md, section 6 of roadmap/README.md) and
# any check that should run what the full run will run. Source it:
#
#   . experiments/wil50_feasibility/fullrun/settings.sh
#   fullrun_check || exit 1
#   latest/bnbdrv -p $FULLRUN_PAR -u $FULLRUN_UB -s $FULLRUN_SYM ... $FULLRUN_INSTANCE
#
# The solver's settings are the [settings] section of fullrun.par (the cpcutt
# arm of pipeline_ab/run_pair.sh, which measured 1.760x in DNN CPU on the
# matched nodes of the node-10 subtree and 2.20x end to end against bestcoret,
# pipeline_ab/runs_cpt/, Nov 7; decided Nov 7: the cut phase in CP form; plus
# the depth spec -x 1:4:1:300:1 of the sampling plan). This file keeps what is
# not a setting: the instance, its symmetry file, the upper bound, the thread
# variables (docs/settings_plan.md, Phase 2, item 7). It exports no policy
# QAP_*/DNN_* variable: the solver rejects the pre-Phase-2 names.

FULLRUN_ROOT=${FULLRUN_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}

# the instance: the transposed wil50 and its symmetry file (sym/README.md); a
# -F i:j spec on wil50.dat is j:i here, and the regional cuts are on side 1
# (cut.side = 1 in fullrun.par)
FULLRUN_INSTANCE=$FULLRUN_ROOT/qaplib/wil50t.dat
FULLRUN_SYM=$FULLRUN_ROOT/latest/SYM/wil50t.sym
FULLRUN_UB=48816
FULLRUN_PAR=$FULLRUN_ROOT/experiments/wil50_feasibility/fullrun/fullrun.par

# one thread per solver process
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
# the one diagnostic the gated arm ran (each assembled pool term checked on 20
# random permutations); a diagnostic, so an environment variable, kept so the
# run is the measured configuration to the letter
export QAP_DNN_POOL_CHECK=20

# true when the files exist and the binary's resolved settings under
# fullrun.par equal the committed golden dump (settings.golden: the
# `settings`/`setting` lines of `bnbdrv --print-settings -p fullrun.par`, less
# the par.file line, whose path depends on the checkout). Any difference fails,
# a changed default included. Regenerate the golden only in a commit that means
# to change the full run's settings:
#   (. fullrun/settings.sh; fullrun_dump latest/bnbdrv > fullrun/settings.golden)
fullrun_dump() {
  "$1" --print-settings -p "$FULLRUN_PAR" 2>/dev/null | grep -E '^settings? ' | grep -v '^setting par\.file '
}
fullrun_check() {
  local b=${1:-$FULLRUN_ROOT/latest/bnbdrv} golden=$FULLRUN_ROOT/experiments/wil50_feasibility/fullrun/settings.golden ok=0
  [ -f "$FULLRUN_INSTANCE" ] || { echo "fullrun: $FULLRUN_INSTANCE is missing (sym/make_sym.py --transpose)"; ok=1; }
  [ -f "$FULLRUN_SYM" ] || { echo "fullrun: $FULLRUN_SYM is missing (sym/README.md)"; ok=1; }
  [ -f "$FULLRUN_PAR" ] || { echo "fullrun: $FULLRUN_PAR is missing"; return 1; }
  [ -x "$b" ] || { echo "fullrun: $b is missing: make solver"; return 1; }
  [ -f "$golden" ] || { echo "fullrun: $golden is missing"; return 1; }
  if ! diff <(fullrun_dump "$b") "$golden" >/dev/null; then
    echo "fullrun: $b's settings under $FULLRUN_PAR differ from $golden:"
    diff <(fullrun_dump "$b") "$golden" | sed 's/^/  /'
    ok=1
  fi
  return $ok
}
