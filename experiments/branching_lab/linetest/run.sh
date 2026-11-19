#!/bin/bash
# The depth-7 line test at tai60b's chain node 44:45,43:58,1:33,6:24,7:19,0:32,10:39 (README.md).
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
BUILD=${BUILD:-/home/user/wt-linetest}
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
  { echo "./bnbdrv $*"; echo "QAP_FORCE_LINE=${QAP_FORCE_LINE:-} QAP_LINE_RANK=${QAP_LINE_RANK:-}";
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
