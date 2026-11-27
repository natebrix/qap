#!/bin/bash
# wil50's current estimate settings (../README.md): the full run's
# (../../wil50_feasibility/fullrun/), its bottom block on rule 9 in place of rule 2
# (experiments/rule9_bottom/). The instance, its symmetry file, the incumbent,
# and the params file whose [settings] section holds every solver setting. Source it:
#   . experiments/estimates/wil50/settings.sh && fullrun_check || exit 1
FULLRUN_ROOT=${FULLRUN_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
FULLRUN_INSTANCE=$FULLRUN_ROOT/qaplib/wil50t.dat
FULLRUN_SYM=$FULLRUN_ROOT/latest/SYM/wil50t.sym
FULLRUN_UB=48816
FULLRUN_PAR=$FULLRUN_ROOT/experiments/estimates/wil50/wil50.par
FULLRUN_GOLDEN=$FULLRUN_ROOT/experiments/estimates/wil50/settings.golden
. "$FULLRUN_ROOT/experiments/estimates/lib.sh"
