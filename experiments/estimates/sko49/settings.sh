#!/bin/bash
# sko49's current estimate settings (../README.md): the instance, its symmetry
# file, the incumbent, and the params file whose [settings] section holds every
# solver setting. Source it:
#   . experiments/estimates/sko49/settings.sh && fullrun_check || exit 1
FULLRUN_ROOT=${FULLRUN_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
FULLRUN_INSTANCE=$FULLRUN_ROOT/qaplib/sko49t.dat
FULLRUN_SYM=$FULLRUN_ROOT/latest/SYM/sko49t.sym
FULLRUN_UB=23386
FULLRUN_PAR=$FULLRUN_ROOT/experiments/estimates/sko49/sko49.par
FULLRUN_GOLDEN=$FULLRUN_ROOT/experiments/estimates/sko49/settings.golden
. "$FULLRUN_ROOT/experiments/estimates/lib.sh"
export QAP_DNN_POOL_CHECK=20   # the full run diagnostic (each assembled pool term checked on 20 random permutations), as in ../../wil50_feasibility/fullrun/settings.sh
