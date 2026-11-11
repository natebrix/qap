#!/bin/bash
# sko42 tuning arm A (README.md here): sko42t with its symmetry, the optimum as the
# incumbent, A.par. Source it, then fullrun_check (../../lib.sh).
FULLRUN_ROOT=${FULLRUN_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)}
FULLRUN_INSTANCE=$FULLRUN_ROOT/qaplib/sko42t.dat
FULLRUN_SYM=$FULLRUN_ROOT/latest/SYM/sko42t.sym
FULLRUN_UB=15812
FULLRUN_PAR=$FULLRUN_ROOT/experiments/estimates/sko42/tune/A.par
FULLRUN_GOLDEN=$FULLRUN_ROOT/experiments/estimates/sko42/tune/A.golden
. "$FULLRUN_ROOT/experiments/estimates/lib.sh"
