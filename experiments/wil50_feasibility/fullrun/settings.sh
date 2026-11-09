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
export DNN_AMX=1 DNN_CHOLQR_PASSES=1 QAP_DNN_CHEB=3 QAP_DNN_TIGHT_LAM=2 QAP_DNN_REFRESH=40
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
export QAP_DNN_SHIP=32 DNN_R32_NONCERT=1
export DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 DNN_ADMM_CROSS=2
export DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1
# one thread per solver process
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

# true when the files exist and the binary knows every switch above
fullrun_check() {
  local b=${1:-$FULLRUN_ROOT/latest/bnbdrv} s ok=0
  [ -f "$FULLRUN_INSTANCE" ] || { echo "fullrun: $FULLRUN_INSTANCE is missing (sym/make_sym.py --transpose)"; ok=1; }
  [ -f "$FULLRUN_SYM" ] || { echo "fullrun: $FULLRUN_SYM is missing (sym/README.md)"; ok=1; }
  [ -x "$b" ] || { echo "fullrun: $b is missing: make solver"; return 1; }
  for s in QAP_DNN_CUT_CP QAP_DNN_CUT_DAMP1 QAP_DNN_SHIP DNN_R32_NONCERT DNN_BETA_SCHED DNN_ADMM_CROSS \
           DNN_CORE_LAZYCERT DNN_CORE_NOSYMV DNN_NOSYMW DNN_CORE_REFRESH_PARTIAL; do
    grep -aq "$s" "$b" || { echo "fullrun: $b does not know $s: git pull, make solver"; ok=1; }
  done
  return $ok
}
