#!/bin/bash
# tho40's current estimate settings (../README.md): the instance, its symmetry
# file, the incumbent, and the params file whose [settings] section holds every
# solver setting. Source it:
#   . experiments/estimates/tho40/settings.sh && fullrun_check || exit 1
FULLRUN_ROOT=${FULLRUN_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
FULLRUN_INSTANCE=$FULLRUN_ROOT/qaplib/tho40t.dat
FULLRUN_SYM=$FULLRUN_ROOT/latest/SYM/tho40t.sym
FULLRUN_UB=240516
FULLRUN_PAR=$FULLRUN_ROOT/experiments/estimates/tho40/tho40.par
FULLRUN_GOLDEN=$FULLRUN_ROOT/experiments/estimates/tho40/settings.golden
. "$FULLRUN_ROOT/experiments/estimates/lib.sh"
