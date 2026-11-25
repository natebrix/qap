#!/bin/bash
# Warm starts and inherited cuts, measured on the children of fixed nodes, for section 5 of
# the DNN+ paper (docs/dnn_plus/root_runs.md, "Warm starts and inherited cuts").
#
#   docs/dnn_plus/warm_pool.sh CASE [CASE ...]        # cases: nug30 tho40 wil50 tai35b (test20: a quick check)
#
# A case is a problem and a parent node (-F fixes; "root" for the root). Each arm bounds the
# parent and then each of its children once (bnbdrv -F ... -d depth+1), one core, with the
# search's own code and the problem's recipe. The parent's solve is the same in every arm of
# an experiment, so every arm sees the same children. The first-solve give-up rule is the
# recipe's dnn.giveup in every arm (no dnn.giveup_warm), so decisions can be compared.
#
#   Experiment A, warm starts (no cut phase anywhere):
#     cold   each child from a cold start                     dnn.warm = 0
#     z      from the parent's multiplier Z                   dnn.warm = 1, dnn.warm_y = 0
#     zy     from Z and the parent's Y, conditioned           dnn.warm = 1, dnn.warm_y = 2
#   Experiment B, inherited cuts (the parent runs a cut phase of ITERS CP-ADMM iterations
#   after its first solve; the children run their first solve only, warm as in zy):
#     nopool the children without the parent's cuts           pool.on = 0
#     pool   the children with the parent's active cuts       pool.on = 1
#
# Environment: ITERS (parent's cut phase, default 2000), ARMS (default "cold z zy nopool pool"),
# OUT (default docs/dnn_plus/runs_warm). Results: OUT/CASE/ARM.{par,log,trace};
# docs/dnn_plus/warm_pool_summary.py OUT prints the tables.
set -uo pipefail
cd "$(dirname "$0")/../.."
R=$(pwd)
ITERS=${ITERS:-2000}
ARMS=${ARMS:-cold z zy nopool pool}
OUT=${OUT:-$R/docs/dnn_plus/runs_warm}
[ -x latest/bnbdrv ] || { echo "build first: make"; exit 2; }
[ $# -ge 1 ] || { sed -n '2,22p' "$0"; exit 2; }
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
REC=infra/phase1/recipes
REG=experiments/bound_lab/regions/regions

for C in "$@"; do
  CUTS=""; GRID=""; FIX=""; ARMS_C=$ARMS
  case $C in
    nug30)  PAR=$REC/nug30.par; DAT=qaplib/nug30.dat; BEST=6124; STEP=2; FIX="0:14"       # depth 1
            GRID="cut.grid = 5:6|cut.side = 1|cut.cp = 100:1000|cut.warmy = 1" ;;         # its recipe has no cut phase
    tho40)  PAR=$REC/tho40c_sb.par; DAT=qaplib/tho40t.dat; BEST=240516; STEP=2; FIX="26:10,8:0,33:32" ;;   # depth 3 (bound_lab/tho40_search, t1's parent)
    wil50)  PAR=$REC/wil50_d2r10.par; DAT=qaplib/wil50t.dat; BEST=48816; STEP=2; FIX="27:39,48:30,49:25" ;;  # depth 3 (pipeline_ab's node 10)
    tai35b) PAR=$REC/tai50b.par; DAT=qaplib/tai35b.dat; BEST=283315445; STEP=1; FIX=""     # the root; experiment B only
            CUTS=$R/$REG/tai35b.comb.txt; ARMS_C=$(echo $ARMS | tr ' ' '\n' | grep -E '^(nopool|pool)$' | tr '\n' ' ') ;;
    test20) PAR=$REC/nug30.par; DAT=qaplib/nug20.dat; BEST=2570; STEP=2; FIX=""           # a quick check of the script
            GRID="cut.grid = 4:5|cut.side = 1|cut.cp = 100:1000|cut.warmy = 1" ;;
    *) echo "unknown case $C"; exit 2 ;;
  esac
  for f in "$PAR" "$DAT" ${CUTS:+"$CUTS"}; do [ -f "$f" ] || { echo "$C: missing $f"; exit 2; }; done
  if [ -z "$FIX" ]; then K=0; else K=$(echo "$FIX" | awk -F, '{print NF}'); fi
  D=$OUT/$C; mkdir -p "$D"
  for ARM in $ARMS_C; do
    grep -v -E '^(dnn\.(depths|warm|warm_y|giveup_warm|retry)|pool\.|sb\.|cut\.(zone|maxdepth|handover|stall|cp_iters|giveup|giveup_depth|freeze|cp_wc_gap|regions|equip|pgate))' \
      "$PAR" > "$D/$ARM.par"
    { echo "dnn.depths = $K:$((K + 1))"
      [ -n "$CUTS" ] && echo "cut.regions = $CUTS"
      [ -n "$GRID" ] && echo "$GRID" | tr '|' '\n'
      case $ARM in
        cold)   echo "dnn.warm = 0"; echo "cut.zone = 0"; echo "pool.on = 0" ;;
        z)      echo "dnn.warm = 1"; echo "dnn.warm_y = 0"; echo "cut.zone = 0"; echo "pool.on = 0" ;;
        zy)     echo "dnn.warm = 1"; echo "dnn.warm_y = 2"; echo "cut.zone = 0"; echo "pool.on = 0" ;;
        nopool|pool)
                echo "dnn.warm = 1"; echo "dnn.warm_y = 2"
                echo "cut.zone = 1"; echo "cut.maxdepth = $K"; echo "cut.stall = 0"; echo "cut.freeze = 0"
                echo "cut.cp_iters = $ITERS"; echo "cut.giveup = 100000000:0"; echo "cut.giveup_depth = 0"
                [ "$ARM" = pool ] && echo "pool.on = 1" || echo "pool.on = 0" ;;
        *) echo "unknown arm $ARM" >&2; exit 2 ;;
      esac; } >> "$D/$ARM.par"
    rm -f "$D/$ARM.trace"
    echo "$(date '+%m-%d %H:%M:%S') $C $ARM: start"
    ( cd latest && QAP_DNN_TRACE="$D/$ARM.trace" \
        ./bnbdrv -p "$D/$ARM.par" -u $((BEST + STEP)) ${FIX:+-F "$FIX"} -d $((K + 1)) "../$DAT" > "$D/$ARM.log" 2>&1 )
    echo "$(date '+%m-%d %H:%M:%S') $C $ARM: exit $? nodes $(wc -l < "$D/$ARM.trace" 2>/dev/null || echo 0)"
  done
  { echo "case $C commit $(git rev-parse --short HEAD) host $(uname -sm | tr ' ' '-') iters $ITERS"; } > "$D/info.txt"
done
