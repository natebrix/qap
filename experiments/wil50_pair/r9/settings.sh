#!/bin/bash
# experiments/wil50_pair, arm r9: the same settings with wil50.par's rule blocks, rule 9 at depth 2 (../README.md).
# The incumbent is the next run's, 48,818. For experiments/estimates/subtree.sh.
FULLRUN_ROOT=${FULLRUN_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}
FULLRUN_INSTANCE=$FULLRUN_ROOT/qaplib/wil50t.dat
FULLRUN_SYM=$FULLRUN_ROOT/latest/SYM/wil50t.sym
FULLRUN_UB=48818
FULLRUN_PAR=$FULLRUN_ROOT/experiments/wil50_pair/r9/r9.par
FULLRUN_GOLDEN=$FULLRUN_ROOT/experiments/wil50_pair/r9/settings.golden
. "$FULLRUN_ROOT/experiments/estimates/lib.sh"
