#!/bin/bash
# lsla (from linetest/run.sh, with a t1 target). The depth-7 line test at tai60b's chain node 44:45,43:58,1:33,6:24,7:19,0:32,10:39 (README.md).
# The design of cognac's experiments/tai60b_lines/run.sh, two levels deeper. Production settings
# (infra/phase1/recipes/tai60b.par, its region file), incumbent 608,215,055 = best known + 1.
#
#   run.sh parent          the depth-6 chain node 44:45,43:58,1:33,6:24,7:19,0:32 cold from -F, branched on
#                          facility 10 as production did (QAP_FORCE_LINE=6:r10); its 54 depth-7 children set
#                          aside in runs/parent/out.state (-d 6: bnbdrv sets aside nodes deeper than -d)
#   run.sh kids            every depth-7 child bounded from the parent's state (-d 7, the params' own rule
#                          where one branches), for the comparison with production's depth-7 outcomes
#   run.sh arm NAME LINE   10:39 warm from the parent's state (cutstate.py), branched on LINE (c54, c47,
#                          r31), its children bounded, the grandchildren set aside (-d 8); checkpoint every
#                          300 CPU s
#   run.sh resume NAME LINE  an arm stopped part-way, from its last checkpoint (traces append)
#   run.sh trim NAME       a stopped arm's checkpoint without the depth-9 nodes and their DNN states
#
# BUILD (default /home/user/wt-linetest) is a worktree of cognac 0e6c2ecf with
# experiments/line_signals/experiment_build.patch applied and `make -C latest` (with scipy-openblas32).
set -euo pipefail
D=$(cd "$(dirname "$0")"; pwd)
BUILD=${BUILD:-/home/user/wt-cognac}
PAR=$BUILD/infra/phase1/recipes/tai60b.par
DAT=$BUILD/qaplib/tai60b.dat
UB=608215055
NODE6=44:45,43:58,1:33,6:24,7:19,0:32
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
grep -q QAP_FORCE_LINE "$BUILD/latest/bnb.C" || { echo "$BUILD/latest is not the experiment build"; exit 2; }
cd "$BUILD/latest"   # cut.regions in tai60b.par is relative to latest/

go() {   # TAG, then bnbdrv's arguments (APPEND=1: add to cmd.txt and the logs, as a resume does)
  local tag=$1; shift
  local out=$D/runs/$tag
  { echo "./bnbdrv $*"; echo "QAP_FORCE_LINE=${QAP_FORCE_LINE:-} QAP_LINE_RANK=${QAP_LINE_RANK:-}"; env | grep '^QAP_LS' || true;
    git -C "$BUILD" rev-parse --short HEAD; date -u +%FT%TZ; } > "$out/cmd.new"
  if [ "${APPEND:-0}" = 1 ]; then cat "$out/cmd.new" >> "$out/cmd.txt"; rm "$out/cmd.new"; else mv "$out/cmd.new" "$out/cmd.txt"; fi
  set +e
  QAP_NODE_TRACE=$out/nodes.txt QAP_DNN_TRACE=$out/dnn.txt ./bnbdrv "$@" >> "$out/log.txt" 2>> "$out/err.txt"
  local rc=$?
  set -e
  { echo "exit $rc"; date -u +%FT%TZ; } >> "$out/cmd.txt"
  return $rc
}

case ${1:-} in
  parent)
    mkdir -p "$D/runs/parent"; rm -f "$D/runs/parent/"{nodes,dnn,log,err}.txt
    QAP_FORCE_LINE=6:r10 QAP_LINE_RANK=6 go parent -p "$PAR" -u $UB -F $NODE6 -d 6 -o "$D/runs/parent/out.state" "$DAT" ;;
  kids)
    [ -s "$D/runs/parent/out.state" ] || { echo "run.sh parent first"; exit 2; }
    mkdir -p "$D/runs/kids"; rm -f "$D/runs/kids/"{nodes,dnn,log,err}.txt
    QAP_LINE_RANK=7 go kids -p "$PAR" -r "$D/runs/parent/out.state" -d 7 -o "$D/runs/kids/out.state" -k 600 "$DAT" ;;
  kidsresume)
    APPEND=1 QAP_LINE_RANK=7 go kids -p "$PAR" -r "$D/runs/kids/out.state" -d 7 -o "$D/runs/kids/out.state" -k 600 "$DAT" ;;
  arm)
    NAME=$2; LINE=$3
    [ -s "$D/runs/parent/out.state" ] || { echo "run.sh parent first"; exit 2; }
    mkdir -p "$D/runs/$NAME"; rm -f "$D/runs/$NAME/"{nodes,dnn,log,err}.txt "$D/runs/$NAME/out.state"
    python3 -I "$D/cutstate.py" "$D/runs/parent/out.state" "$D/runs/$NAME/in.state" 7 10:39
    QAP_FORCE_LINE=7:$LINE QAP_LINE_RANK=7 go "$NAME" -p "$PAR" -r "$D/runs/$NAME/in.state" -d 8 \
      -o "$D/runs/$NAME/out.state" -k 300 "$DAT"
    rm -f "$D/runs/$NAME/in.state" ;;
  strong)
    # as arm, under strong.par: tai60b.par with cut.giveup_depth = 8 and cut.cp_iters = 8000, so 10:39's own
    # cut phase runs up to 8,000 iterations instead of stopping at 1,000 (production's settings below it)
    NAME=s$2; LINE=$2
    mkdir -p "$D/runs/$NAME"; rm -f "$D/runs/$NAME/"{nodes,dnn,log,err}.txt "$D/runs/$NAME/out.state"
    python3 -I "$D/cutstate.py" "$D/runs/parent/out.state" "$D/runs/$NAME/in.state" 7 10:39
    QAP_FORCE_LINE=7:$LINE QAP_LINE_RANK=7 go "$NAME" -p "$D/strong.par" --allow-settings-change -r "$D/runs/$NAME/in.state" -d 8 \
      -o "$D/runs/$NAME/out.state" -k 300 "$DAT"
    rm -f "$D/runs/$NAME/in.state" ;;
  t1)
    # lsla T1: 10:39 warm from the parent's state under strong.par (as linetest's strong arms), then QAP_LS mode b
    # at depth 7: NAME, CAPS (QAP_LS_CAP), K; the children set aside (-d 7), not bounded
    NAME=$2; CAPS=$3; K=$4
    mkdir -p "$D/runs/$NAME"; rm -f "$D/runs/$NAME/"{nodes,dnn,log,err,ls,dump}.txt "$D/runs/$NAME/out.state"
    python3 -I "$D/cutstate.py" "$D/runs/parent/out.state" "$D/runs/$NAME/in.state" 7 10:39
    QAP_LS=b:$K:7:7 QAP_LS_CAP=$CAPS QAP_LS_ADD=c47,c54 QAP_LS_SEARCH=${SEARCH:-30:3000:2.5:300} \
      QAP_LS_PERM=$BUILD/QAPSOLN/tai60b.sln QAP_LS_LOG=$D/runs/$NAME/ls.txt QAP_LS_DUMP=$D/runs/$NAME/dump.txt \
      go "$NAME" -p "$D/strong.par" --allow-settings-change -r "$D/runs/$NAME/in.state" -d 7 "$DAT"
    rm -f "$D/runs/$NAME/in.state" ;;
  set8)
    # 10:39 warm from the parent's state under strong.par, branched on LINE, its depth-8 children set aside
    # unbounded (-d 7), for bounding one depth-8 node under different lines (run.sh at8)
    LINE=$2; NAME=set8_$LINE
    mkdir -p "$D/runs/$NAME"; rm -f "$D/runs/$NAME/"{nodes,dnn,log,err}.txt "$D/runs/$NAME/out.state"
    python3 -I "$D/cutstate.py" "$D/runs/parent/out.state" "$D/runs/$NAME/in.state" 7 10:39
    QAP_FORCE_LINE=7:$LINE QAP_LINE_RANK=7 go "$NAME" -p "$D/strong.par" --allow-settings-change -r "$D/runs/$NAME/in.state" -d 7 \
      -o "$D/runs/$NAME/out.state" "$DAT"
    rm -f "$D/runs/$NAME/in.state" ;;
  at8)
    # one depth-8 node FIX (e.g. 14:54) from runs/set8_LINE7/out.state, bounded under strong.par and branched on
    # LINE8 (QAP_FORCE_LINE=8:LINE8), its depth-9 children bounded, depth 10 set aside
    LINE7=$2; FIX=$3; LINE8=$4; NAME=at8_${FIX//:/_}_$LINE8
    mkdir -p "$D/runs/$NAME"; rm -f "$D/runs/$NAME/"{nodes,dnn,log,err}.txt "$D/runs/$NAME/out.state"
    python3 -I "$D/cutstate.py" "$D/runs/set8_$LINE7/out.state" "$D/runs/$NAME/in.state" 8 "$FIX"
    QAP_FORCE_LINE=8:$LINE8 QAP_LINE_RANK=8 go "$NAME" -p "$D/strong.par" --allow-settings-change -r "$D/runs/$NAME/in.state" -d 9 \
      -o "$D/runs/$NAME/out.state" -k 600 "$DAT"
    rm -f "$D/runs/$NAME/in.state" ;;
  set9)
    # 14:54 from runs/set8_c54/out.state under strong.par, branched on LINE8 (default c2, production's), its
    # depth-9 children set aside unbounded (-d 8), for bounding the chain node 24:2 under different lines (at9)
    LINE8=${2:-c2}; NAME=set9_$LINE8
    mkdir -p "$D/runs/$NAME"; rm -f "$D/runs/$NAME/"{nodes,dnn,log,err}.txt "$D/runs/$NAME/out.state"
    python3 -I "$D/cutstate.py" "$D/runs/set8_c54/out.state" "$D/runs/$NAME/in.state" 8 14:54
    QAP_FORCE_LINE=8:$LINE8 QAP_LINE_RANK=8 go "$NAME" -p "$D/strong.par" --allow-settings-change -r "$D/runs/$NAME/in.state" -d 8 \
      -o "$D/runs/$NAME/out.state" "$DAT"
    rm -f "$D/runs/$NAME/in.state" ;;
  at9)
    # LINE8 FIX LINE9: one depth-9 node FIX (e.g. 24:2) from runs/set9_LINE8/out.state, strong.par, branched on
    # LINE9, its depth-10 children bounded, depth 11 set aside
    LINE8=$2; FIX=$3; LINE9=$4; NAME=at9_${FIX//:/_}_$LINE9
    mkdir -p "$D/runs/$NAME"; rm -f "$D/runs/$NAME/"{nodes,dnn,log,err}.txt "$D/runs/$NAME/out.state"
    python3 -I "$D/cutstate.py" "$D/runs/set9_$LINE8/out.state" "$D/runs/$NAME/in.state" 9 "$FIX"
    QAP_FORCE_LINE=9:$LINE9 QAP_LINE_RANK=9 go "$NAME" -p "$D/strong.par" --allow-settings-change -r "$D/runs/$NAME/in.state" -d 10 \
      -o "$D/runs/$NAME/out.state" -k 600 "$DAT"
    rm -f "$D/runs/$NAME/in.state" ;;
  t1d9)
    # mode b at 24:2 (depth 9) from runs/set9_c2/out.state, strong.par: NAME CAPS K; ADD r8 (rule 2's line);
    # rule 9's r3 joins as the rule's own line; children set aside (-d 9), not bounded
    NAME=$2; CAPS=$3; K=$4
    mkdir -p "$D/runs/$NAME"; rm -f "$D/runs/$NAME/"{nodes,dnn,log,err,ls,dump}.txt "$D/runs/$NAME/out.state"
    python3 -I "$D/cutstate.py" "$D/runs/set9_c2/out.state" "$D/runs/$NAME/in.state" 9 24:2
    QAP_LS=b:$K:9:9 QAP_LS_CAP=$CAPS QAP_LS_ADD=${ADD:-r8} QAP_LS_SEARCH=${SEARCH:-30:3000:2.5:300} \
      QAP_LS_PERM=$BUILD/QAPSOLN/tai60b.sln QAP_LS_LOG=$D/runs/$NAME/ls.txt QAP_LS_DUMP=$D/runs/$NAME/dump.txt \
      go "$NAME" -p "$D/strong.par" --allow-settings-change -r "$D/runs/$NAME/in.state" -d 9 "$DAT"
    rm -f "$D/runs/$NAME/in.state" ;;
  root50la)
    # mode b at tai50b's root (production recipe, rule 9), children set aside (-d 0): NAME CAPS K; ADD lines for the
    # lookahead to judge besides the signal's top K and rule 9's two lines; SEARCH as QAP_LS_SEARCH
    NAME=$2; CAPS=$3; K=$4
    mkdir -p "$D/runs/$NAME"; rm -f "$D/runs/$NAME/"{nodes,dnn,log,err,ls,dump}.txt "$D/runs/$NAME/out.state"
    QAP_LS=b:$K:0:0 QAP_LS_CAP=$CAPS QAP_LS_ADD=${ADD:-} QAP_LS_SEARCH=${SEARCH:-30:0:2.5:300:100} \
      QAP_LS_PERM=$BUILD/QAPSOLN/tai50b.sln QAP_LS_LOG=$D/runs/$NAME/ls.txt QAP_LS_DUMP=$D/runs/$NAME/dump.txt \
      go "$NAME" -p "$BUILD/infra/phase1/recipes/tai50b.par" -u 458821518 -d 0 "$BUILD/qaplib/tai50b.dat" ;;
  root50lines)
    # plain QAP_LA (the same per-child lookahead as mode b) at tai50b's root over exactly LINES (QAP_LA_LINES),
    # CAP iterations a child, children set aside (-d 0): NAME CAP LINES. Splits root50la's shortlist across
    # processes, so a container restart costs less; each part carries r31 as a cross-check.
    NAME=$2; CAP=$3; LINES=$4
    mkdir -p "$D/runs/$NAME"; rm -f "$D/runs/$NAME/"{nodes,dnn,log,err,la,children}.txt
    QAP_LA=1:$CAP:0:0 QAP_LA_LINES=$LINES QAP_LA_LOG=$D/runs/$NAME/la.txt QAP_LA_PROGRESS=$D/runs/$NAME/children.txt \
      go "$NAME" -p "$BUILD/infra/phase1/recipes/tai50b.par" -u 458821518 -d 0 "$BUILD/qaplib/tai50b.dat" ;;
  root50)
    # tai50b's root under its production recipe (infra/phase1/recipes/tai50b.par, rule 9), branched on LINE,
    # its depth-1 children bounded, depth 2 set aside; incumbent best known + 1, as the laptop's run
    LINE=$2; NAME=root50_$LINE
    mkdir -p "$D/runs/$NAME"; rm -f "$D/runs/$NAME/"{nodes,dnn,log,err}.txt "$D/runs/$NAME/out.state"
    QAP_FORCE_LINE=0:$LINE QAP_LINE_RANK=0 go "$NAME" -p "$BUILD/infra/phase1/recipes/tai50b.par" -u 458821518 -d 1 \
      -o "$D/runs/$NAME/out.state" "$BUILD/qaplib/tai50b.dat" ;;
  deep)
    # one level deeper, as linetest's deep target: the depth-9 nodes under a strong arm's survivors FIXES (e.g.
    # 11:59 or 13:59,3:59), from runs/s<LINE>/out.state, production's settings and line rule, depth 10 set aside
    LINE=$2; NAME=d9${LINE}_${3//[:,]/_}
    mkdir -p "$D/runs/$NAME"; rm -f "$D/runs/$NAME/"{nodes,dnn,log,err}.txt "$D/runs/$NAME/out.state"
    python3 -I "$D/splitstate.py" "$D/runs/s$LINE/out.state" "$D/runs/$NAME/in.state" 9 "$3"
    go "$NAME" -p "$PAR" --allow-settings-change -r "$D/runs/$NAME/in.state" -d 9 -o "$D/runs/$NAME/out.state" -k 600 "$DAT"
    rm -f "$D/runs/$NAME/in.state" ;;
  sresume)
    NAME=s$2; LINE=$2
    grep -q '^exit 0$' "$D/runs/$NAME/cmd.txt" && { echo "$NAME has finished"; exit 2; }
    QAP_FORCE_LINE=7:$LINE QAP_LINE_RANK=7 APPEND=1 go "$NAME" -p "$D/strong.par" --allow-settings-change -r "$D/runs/$NAME/out.state" -d 8 \
      -o "$D/runs/$NAME/out.state" -k 300 "$DAT" ;;
  resume)
    NAME=$2; LINE=$3
    grep -q '^exit 0$' "$D/runs/$NAME/cmd.txt" && { echo "$NAME has finished"; exit 2; }
    [ -s "$D/runs/$NAME/out.state" ] || { echo "$NAME has no checkpoint: run.sh arm $NAME $LINE"; exit 2; }
    QAP_FORCE_LINE=7:$LINE QAP_LINE_RANK=7 APPEND=1 go "$NAME" -p "$PAR" -r "$D/runs/$NAME/out.state" -d 8 \
      -o "$D/runs/$NAME/out.state" -k 300 "$DAT" ;;
  trim)
    NAME=$2; S=$D/runs/$NAME/out.state
    for p in $(pgrep -x bnbdrv); do
      tr '\0' ' ' < /proc/$p/cmdline | grep -q "runs/$NAME/out.state" && { echo "$NAME is running"; exit 2; }
    done
    k=1; while [ -e "$D/runs/$NAME/setaside.$k.state" ]; do k=$((k+1)); done
    mv "$S" "$D/runs/$NAME/setaside.$k.state"
    python3 -I "$D/trimstate.py" "$D/runs/$NAME/setaside.$k.state" "$S" 8 ;;
  *) awk 'NR > 1 && /^set /{exit} NR > 1' "$0"; exit 2 ;;
esac
