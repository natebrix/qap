#!/bin/bash
# sko42's current estimate settings (../README.md): the instance, its symmetry
# file, the incumbent, and the params file whose [settings] section holds every
# solver setting. Source it:
#   . experiments/estimates/sko42/settings.sh && fullrun_check || exit 1
FULLRUN_ROOT=${FULLRUN_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
FULLRUN_INSTANCE=$FULLRUN_ROOT/qaplib/sko42t.dat
FULLRUN_SYM=$FULLRUN_ROOT/latest/SYM/sko42t.sym
FULLRUN_UB=15812
FULLRUN_PAR=$FULLRUN_ROOT/experiments/estimates/sko42/sko42.par
FULLRUN_GOLDEN=$FULLRUN_ROOT/experiments/estimates/sko42/settings.golden
. "$FULLRUN_ROOT/experiments/estimates/lib.sh"
