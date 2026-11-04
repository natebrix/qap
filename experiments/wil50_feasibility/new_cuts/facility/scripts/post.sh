#!/bin/bash
# All facility-side diagnostics on one dumped iterate: post.sh SPEC Y.bin TAG [PRIMAL]
cd "$(dirname "$0")"
S=$1; Y=$2; T=$3; P=${4:+--primal $4}
L=../logs
taskset -c 3 python3 ydiag.py $S $Y $P > $L/ydiag_$T.log 2>&1
taskset -c 3 python3 box_sep.py $S $Y --starts 60 > $L/box_$T.log 2>&1
taskset -c 3 python3 facreg_search.py $S $Y --triples 120 --starts 6 > $L/facreg_$T.log 2>&1
taskset -c 3 python3 sep_dist.py $S $Y --fam hull3 --top 300 2>&1 | grep -v Warn > $L/hull3_$T.log
taskset -c 3 python3 disp_check.py $S $Y > $L/disp_$T.log 2>&1
echo "done $T"
