#!/bin/bash
# Root bounds for the DNN+ paper (docs/dnn_plus/root_runs.md): for each PROBLEM, the certified
# DNN bound without cuts ("plain") and with regional cuts ("plus"), root only, one core each.
#
#   docs/dnn_plus/root_bound.sh PROBLEM [PROBLEM ...]
#   printf '%s\n' sko42 sko49 wil50 tho40 | xargs -P 4 -n 1 docs/dnn_plus/root_bound.sh   # four at a time
#
# Problems: sko42 sko49 wil50 tho40 nug30 tai25b tai30b tai35b tai40b tai50b tai60b.
# Environment: ITERS (cut-phase iteration limit, default 6000), ARMS ("plain plus"),
# OUT (default docs/dnn_plus/runs).
#
# Each arm is the search's own code (latest/bnbdrv -d 0) under the problem's recipe, with these
# changes, written to OUT/PROBLEM/ARM.par: the DNN bound at depth 0 only; no strong branching;
#   plain: no cut phase, no give-up: the first solve runs to its limit (4,000 iterations);
#   plus:  the recipe's first solve, then the cut phase at any gap (cut.zone = 1) as one
#          CP-ADMM solve of ITERS iterations, with no stall, hand-over or give-up rule.
# The incumbent is the best known value plus one objective step, so nothing fathoms.
# Results: a line per arm in OUT/results.txt; OUT/PROBLEM/ARM.{log,trace}; the plus trace carries
# the certified value at every check (| hist it:value ...), for the stopping rule and the figure.
set -uo pipefail
cd "$(dirname "$0")/../.."
R=$(pwd)
ITERS=${ITERS:-6000}
ARMS=${ARMS:-plain plus}
OUT=${OUT:-$R/docs/dnn_plus/runs}
[ -x latest/bnbdrv ] || { echo "build first: make"; exit 2; }
[ $# -ge 1 ] || { sed -n '2,10p' "$0"; exit 2; }
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
REC=infra/phase1/recipes
REG=experiments/bound_lab/regions/regions

for P in "$@"; do
  CUTS=""; GRID=""
  case $P in
    sko42)  PAR=$REC/sko42_d2r10_sb.par; DAT=qaplib/sko42t.dat; BEST=15812;     STEP=2 ;;
    sko49)  PAR=$REC/sko49_d2r10_sb.par; DAT=qaplib/sko49t.dat; BEST=23386;     STEP=2 ;;
    wil50)  PAR=$REC/wil50_d2r10.par;    DAT=qaplib/wil50t.dat; BEST=48816;     STEP=2 ;;
    tho40)  PAR=$REC/tho40c_sb.par;      DAT=qaplib/tho40t.dat; BEST=240516;    STEP=2 ;;
    nug30)  PAR=$REC/nug30.par;          DAT=qaplib/nug30.dat;  BEST=6124;      STEP=2
            GRID="pool.on = 1|cut.grid = 5:6|cut.side = 1|cut.cp = 100:1000|cut.warmy = 1" ;;   # its recipe has no cut phase
    tai25b|tai30b|tai35b|tai40b)
            PAR=$REC/tai50b.par; DAT=qaplib/$P.dat; STEP=1; CUTS=$R/$REG/$P.comb.txt
            BEST=$(awk 'NR == 1 {print $2; exit}' QAPSOLN/$P.sln | tr -d '\r[:space:]') ;;
    tai50b) PAR=$REC/tai50b.par; DAT=qaplib/tai50b.dat; BEST=458821517; STEP=1; CUTS=$R/$REC/tai50b.regions ;;
    tai60b) PAR=$REC/tai60b.par; DAT=qaplib/tai60b.dat; BEST=608215054; STEP=1; CUTS=$R/$REC/tai60b.regions ;;
    *) echo "unknown problem $P"; exit 2 ;;
  esac
  for f in "$PAR" "$DAT" ${CUTS:+"$CUTS"}; do [ -f "$f" ] || { echo "$P: missing $f"; exit 2; }; done
  [[ $BEST =~ ^[0-9]+$ ]] || { echo "$P: no best known value"; exit 2; }
  D=$OUT/$P; mkdir -p "$D"
  for ARM in $ARMS; do
    # the recipe without the keys this script sets
    grep -v -E '^(dnn\.depths|sb\.|cut\.(zone|handover|stall|cp_iters|giveup|giveup_depth|freeze|cp_wc_gap|regions|equip|pgate))' \
      "$PAR" > "$D/$ARM.par"
    { echo "dnn.depths = 0:0"
      [ -n "$CUTS" ] && echo "cut.regions = $CUTS"
      [ -n "$GRID" ] && [ "$ARM" = plus ] && echo "$GRID" | tr '|' '\n'
      if [ "$ARM" = plain ]; then
        echo "cut.zone = 0"; grep -v '^dnn\.giveup' "$D/$ARM.par" > "$D/$ARM.tmp"; mv "$D/$ARM.tmp" "$D/$ARM.par"
        echo "dnn.giveup = 0"
      else
        echo "cut.zone = 1"; echo "cut.stall = 0"; echo "cut.freeze = 0"; echo "cut.cp_iters = $ITERS"
        echo "cut.giveup = 100000000:0"; echo "cut.giveup_depth = 0"    # never stops; makes the phase record its history
      fi; } >> "$D/$ARM.par"
    rm -f "$D/$ARM.trace"
    echo "$(date '+%m-%d %H:%M:%S') $P $ARM: start"
    ( cd latest && QAP_DNN_TRACE="$D/$ARM.trace" QAP_DNN_TRACE_HIST=1 \
        ./bnbdrv -p "$D/$ARM.par" -u $((BEST + STEP)) -d 0 "../$DAT" > "$D/$ARM.log" 2>&1 )
    rc=$?
    [ -s "$D/$ARM.trace" ] || { echo "$P $ARM: no trace (exit $rc); see $D/$ARM.log"; continue; }
    awk -v p="$P" -v arm="$ARM" -v best="$BEST" -v c="$(git rev-parse --short HEAD)" -v h="$(uname -sm | tr ' ' '-')" '
      { for (i = 1; i <= NF; i++) { if ($i == "lb" && $(i-1) == "first") lb = $(i+1); if ($i == "best") cb = $(i+1)
          if ($i == "iters") { n++; it[n] = $(i+1) } if ($i == "cpu") { m++; cpu[m] = $(i+1) } } }
      END { v = (arm == "plus" && cb != "") ? cb : lb; its = it[1] + (arm == "plus" ? it[2] : 0); t = cpu[1] + (arm == "plus" ? cpu[2] : 0)
            printf "%s %s bound %.4f gap %.4f%% iters %d cpu_s %.0f first_lb %.4f commit %s host %s\n", p, arm, v, 100 * (best - v) / best, its, t, lb, c, h }' \
      "$D/$ARM.trace" | tee -a "$OUT/results.txt"
  done
done
