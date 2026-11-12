#!/bin/bash
# sko49 arm Cw, C with cut.zone 0.018 (README.md here): sko49t with its symmetry, the best known value as the
# incumbent, Cw.par. Source it, then fullrun_check (../../lib.sh).
FULLRUN_ROOT=${FULLRUN_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)}
FULLRUN_INSTANCE=$FULLRUN_ROOT/qaplib/sko49t.dat
FULLRUN_SYM=$FULLRUN_ROOT/latest/SYM/sko49t.sym
FULLRUN_UB=23386
FULLRUN_PAR=$FULLRUN_ROOT/experiments/estimates/sko49/tune/Cw.par
FULLRUN_GOLDEN=$FULLRUN_ROOT/experiments/estimates/sko49/tune/Cw.golden
. "$FULLRUN_ROOT/experiments/estimates/lib.sh"
