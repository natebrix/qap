#!/bin/bash
# The full wil50 run's settings: one place for the launcher (not built yet),
# the sampling plan (roadmap/tree_size.md, section 6 of roadmap/README.md) and
# any check that should run what the full run will run. Source it:
#
#   . experiments/wil50_feasibility/fullrun/settings.sh
#   fullrun_check || exit 1
#
# It exports the solver environment and sets FULLRUN_INSTANCE, FULLRUN_SYM,
# FULLRUN_UB and FULLRUN_ROOT. The command-line depth spec (-x) and the
# frontier are the launcher's (README.md here). This is the `cpcutt` arm of
# pipeline_ab/run_pair.sh, which measured 1.760x in DNN CPU on the matched
# nodes of the node-10 subtree and 2.20x end to end against `bestcoret`
# (pipeline_ab/runs_cpt/, Nov 7). Decided Nov 7: the cut phase in CP form.

FULLRUN_ROOT=${FULLRUN_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}

# the instance: the transposed wil50 and its symmetry file (sym/README.md); a
# -F i:j spec on wil50.dat is j:i here, and the regional cuts are on side 1
FULLRUN_INSTANCE=$FULLRUN_ROOT/qaplib/wil50t.dat
FULLRUN_SYM=$FULLRUN_ROOT/latest/SYM/wil50t.sym
FULLRUN_UB=48816

# the DNN bound (the base of every pipeline_ab arm)
export QAP_DNN_GAMMA1=1.0 QAP_DNN_GAMMA=1.0 QAP_DNN_EXT=20
export DNN_AMX=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40
export QAP_DNN_WARM_Y=2 QAP_DNN_POOL=1 QAP_DNN_POOL_CHECK=20
# Pinned Nov 9: 34e347c (Nov 8) changed these two defaults (ADMM beta 0.3 ->
# 0.21, give-up confirmation zone 0 -> 0.003). The gated cpcutt arm and the
# sampling plan's stages A and B ran the old values with neither variable set.
# Adopting the new values for wil50 is an A/B, not a default change
# (docs/settings_plan.md, step 0).
export QAP_DNN_BETA=0.3 QAP_DNN_GU_ZONE=0
# the cut phase: regional triangles on the 5:10 grid of the second matrix,
# in CP form (QAP_DNN_CUT_CP: separate 1000 every 100 iterations inside one
# CP-ADMM solve; roadmap/admm/cpadmm/README.md, "In the search" and "Side 1")
export QAP_DNN_CUT_GRID=5:10 QAP_DNN_CUT_SIDE=1
export QAP_DNN_CUT_ZONE=0.005 QAP_DNN_CUT_WARMY=1 QAP_DNN_CUT_ROUNDS=12 QAP_DNN_CUT_ITERS=200
export QAP_DNN_CUT_FREEZE=1000 QAP_DNN_CUT_DAMP1=0.3
export QAP_DNN_CUT_CP=100:1000
# the parent's state in state files, fp32 for the non-certificate iterations,
# the warm-call beta schedule, the crossing certificate, the core switches
export QAP_DNN_SHIP=32
export DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 DNN_ADMM_CROSS=2
export DNN_CORE_LAZYCERT=1
# DNN_CHOLQR_PASSES=1, DNN_R32_NONCERT=1, DNN_CORE_NOSYMV=1, DNN_NOSYMW=1 and
# DNN_CORE_REFRESH_PARTIAL=1 are the solver's only behaviour since Nov 9 (the
# variables were removed with the losing paths, docs/settings_plan.md step 1).
# one thread per solver process
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

# true when the files exist and the binary's resolved settings under this
# environment equal the committed golden dump (settings.golden: the
# `settings`/`setting` lines of `bnbdrv --print-settings`, docs/settings_plan.md
# Phase 1). Any difference fails, a changed default included. Regenerate the
# golden only in a commit that means to change the full run's settings:
#   (. fullrun/settings.sh; latest/bnbdrv --print-settings | grep -E '^settings? ' > fullrun/settings.golden)
# The launcher's own -p and -x are not in the golden (they show as defaults);
# a launcher that passes them should diff its own dump the same way.
fullrun_check() {
  local b=${1:-$FULLRUN_ROOT/latest/bnbdrv} golden=$FULLRUN_ROOT/experiments/wil50_feasibility/fullrun/settings.golden ok=0
  [ -f "$FULLRUN_INSTANCE" ] || { echo "fullrun: $FULLRUN_INSTANCE is missing (sym/make_sym.py --transpose)"; ok=1; }
  [ -f "$FULLRUN_SYM" ] || { echo "fullrun: $FULLRUN_SYM is missing (sym/README.md)"; ok=1; }
  [ -x "$b" ] || { echo "fullrun: $b is missing: make solver"; return 1; }
  [ -f "$golden" ] || { echo "fullrun: $golden is missing"; return 1; }
  if ! diff <("$b" --print-settings 2>/dev/null | grep -E '^settings? ') "$golden" >/dev/null; then
    echo "fullrun: $b's settings differ from $golden:"
    diff <("$b" --print-settings 2>/dev/null | grep -E '^settings? ') "$golden" | sed 's/^/  /'
    ok=1
  fi
  return $ok
}
