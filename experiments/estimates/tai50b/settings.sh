#!/bin/bash
# tai50b's current estimate settings (../README.md): the instance, its symmetry
# file, the incumbent, and the params file whose [settings] section holds every
# solver setting. Source it:
#   . experiments/estimates/tai50b/settings.sh && fullrun_check || exit 1
FULLRUN_ROOT=${FULLRUN_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
FULLRUN_INSTANCE=$FULLRUN_ROOT/qaplib/tai50b.dat
FULLRUN_SYM=
FULLRUN_UB=458821517
FULLRUN_PAR=$FULLRUN_ROOT/experiments/estimates/tai50b/tai50b.par
FULLRUN_GOLDEN=$FULLRUN_ROOT/experiments/estimates/tai50b/settings.golden
. "$FULLRUN_ROOT/experiments/estimates/lib.sh"
