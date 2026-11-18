#!/bin/bash
# Dump what the branching rules read at the line-truth nodes (README.md): QPB's U, the bound, the cutoff,
# the DNN marginals and every line's lookahead score, with the experiment build's QAP_LINE_SCORES
# (line_signals.patch on the tai60b line study's build; never cognac's latest/).
#   dump.sh tho30 NODE | tho40 NODE | tai60b | tai60b_la
set -euo pipefail
D=$(cd "$(dirname "$0")"; pwd)
R=$(cd "$D/../.."; pwd)
BUILD=${BUILD:-/home/user/wt-lines}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
grep -q QAP_LINE_SCORES "$BUILD/latest/branch.C" || { echo "$BUILD/latest is not the experiment build"; exit 2; }
cd "$BUILD/latest"
mkdir -p "$D/dumps"
case $1 in
  tho30)
    tag=tho30_${2//[:,]/-}; depth=$(($(tr -cd ':' <<< "$2" | wc -c)))
    QAP_LINE_SCORES=$depth QAP_LINE_SCORES_LA=1 QAP_XBAR_DEBUG=2 ./bnbdrv \
      -p "$R/experiments/tho30_solve/line_truth/cutd2_cp_r9.par" -u 149936 -F "$2" -d $depth -x 1:7:1:200:1 \
      "$R/qaplib/tho30.dat" > "$D/dumps/$tag.out" 2> "$D/dumps/$tag.err" ;;
  tho40)
    tag=tho40_${2//[:,]/-}; depth=$(($(tr -cd ':' <<< "$2" | wc -c)))
    QAP_LINE_SCORES=$depth QAP_LINE_SCORES_LA=1 QAP_XBAR_DEBUG=2 ./bnbdrv \
      -p "$R/infra/phase1/recipes/nug30.par" -u 240516 -F "$2" -d $depth "$R/qaplib/tho40t.dat" \
      > "$D/dumps/$tag.out" 2> "$D/dumps/$tag.err" ;;
  tai60b)   # 43:58 warm from 44:45's state, as the line study's arms (../tai60b_lines/run.sh arm)
    L=$R/experiments/tai60b_lines
    python3 -I "$R/experiments/bound_lab/tho40_search/cutstate.py" "$L/runs/parent/out.state" "$D/dumps/tai60b.in.state" 2 43:58
    QAP_LINE_SCORES=2 QAP_LINE_SCORES_LA=1 QAP_XBAR_DEBUG=2 ./bnbdrv -p "$R/infra/phase1/recipes/tai60b.par" \
      -r "$D/dumps/tai60b.in.state" -d 2 "$R/qaplib/tai60b.dat" > "$D/dumps/tai60b_43-58.out" 2> "$D/dumps/tai60b_43-58.err"
    rm -f "$D/dumps/tai60b.in.state" ;;
  tai60b_la)   # the same 43:58, branched by the DNN lookahead over the line study's four lines (README.md, 5)
    L=$R/experiments/tai60b_lines
    python3 -I "$R/experiments/bound_lab/tho40_search/cutstate.py" "$L/runs/parent/out.state" "$D/dumps/tai60b_la.in.state" 2 43:58
    QAP_LA=1:${CAP:-30}:2:2 QAP_LA_LINES=r1,r45,c26,r8 QAP_LA_LOG=$D/dumps/tai60b_la.log QAP_XBAR_DEBUG=1 ./bnbdrv \
      -p "$R/infra/phase1/recipes/tai60b.par" -r "$D/dumps/tai60b_la.in.state" -d 2 "$R/qaplib/tai60b.dat" \
      > "$D/dumps/tai60b_la.out" 2> "$D/dumps/tai60b_la.err"
    rm -f "$D/dumps/tai60b_la.in.state" ;;
esac
