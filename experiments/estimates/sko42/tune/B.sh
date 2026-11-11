#!/bin/bash
# sko42 tuning arm B (README.md here): sko42t with its symmetry, the optimum as the
# incumbent, B.par. Source it, then fullrun_check (../../lib.sh).
FULLRUN_ROOT=${FULLRUN_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)}
FULLRUN_INSTANCE=$FULLRUN_ROOT/qaplib/sko42t.dat
FULLRUN_SYM=$FULLRUN_ROOT/latest/SYM/sko42t.sym
FULLRUN_UB=15812
FULLRUN_PAR=$FULLRUN_ROOT/experiments/estimates/sko42/tune/B.par
FULLRUN_GOLDEN=$FULLRUN_ROOT/experiments/estimates/sko42/tune/B.golden
. "$FULLRUN_ROOT/experiments/estimates/lib.sh"
