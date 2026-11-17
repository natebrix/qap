#!/bin/bash
# The tai60b line study at 43:58 (README.md): which line should the depth-2 node on the optimum's path
# (44:45, 43:58) branch on? Production settings (infra/phase1/recipes/tai60b.par), incumbent 608,215,055.
#
#   run.sh parent          44:45 cold from -F, branched on facility 43 as production did (QAP_FORCE_LINE
#                          1:r43); its 59 children set aside in runs/parent/out.state (-d 1 -o)
#   run.sh arm NAME LINE   43:58 warm from 44:45's state (cutstate.py), branched on LINE (r1, c26, ...;
#                          'rule9' for the params' own choice), its children bounded warm from it, the
#                          grandchildren of those that do not fathom set aside (-d 3 -o); a checkpoint
#                          in runs/NAME/out.state every 900 CPU s (-k)
#   run.sh resume NAME LINE  an arm killed part-way (no "exit" in its cmd.txt) from its last checkpoint;
#                          its traces go on in the same files, so the nodes after that checkpoint come twice
#                          (score.py keeps a node's last line)
#
# Needs the experiment build: BUILD (default /home/user/wt-lines) is a worktree of cognac with
# force_line.patch applied and `make -C latest`; production's latest/ does not have QAP_FORCE_LINE.
set -euo pipefail
D=$(cd "$(dirname "$0")"; pwd)
R=$(cd "$D/../.."; pwd)
BUILD=${BUILD:-/home/user/wt-lines}
PAR=$R/infra/phase1/recipes/tai60b.par
UB=608215055
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
    mkdir -p "$D/runs/parent"; rm -f "$D/runs/parent/"{nodes,dnn}.txt
    QAP_FORCE_LINE=1:r43 QAP_LINE_RANK=1 go parent -p "$PAR" -u $UB -F 44:45 -d 1 -o "$D/runs/parent/out.state" \
      "$R/qaplib/tai60b.dat" ;;
  arm)
    NAME=$2; LINE=$3
    [ -s "$D/runs/parent/out.state" ] || { echo "run.sh parent first"; exit 2; }
    mkdir -p "$D/runs/$NAME"; rm -f "$D/runs/$NAME/"{nodes,dnn,log,err}.txt "$D/runs/$NAME/out.state"
    python3 -I "$R/experiments/bound_lab/tho40_search/cutstate.py" "$D/runs/parent/out.state" "$D/runs/$NAME/in.state" 2 43:58
    F=""; [ "$LINE" = rule9 ] || F=2:$LINE
    QAP_FORCE_LINE=$F QAP_LINE_RANK=2 go "$NAME" -p "$PAR" -r "$D/runs/$NAME/in.state" -d 3 \
      -o "$D/runs/$NAME/out.state" -k 900 "$R/qaplib/tai60b.dat"
    rm -f "$D/runs/$NAME/in.state" ;;
  resume)
    NAME=$2; LINE=$3
    grep -q '^exit' "$D/runs/$NAME/cmd.txt" && { echo "$NAME has finished"; exit 2; }
    [ -s "$D/runs/$NAME/out.state" ] || { echo "$NAME has no checkpoint: run.sh arm $NAME $LINE"; exit 2; }
    F=""; [ "$LINE" = rule9 ] || F=2:$LINE
    QAP_FORCE_LINE=$F QAP_LINE_RANK=2 APPEND=1 go "$NAME" -p "$PAR" -r "$D/runs/$NAME/out.state" -d 3 \
      -o "$D/runs/$NAME/out.state" -k 900 "$R/qaplib/tai60b.dat"
    rm -f "$D/runs/$NAME/in.state" ;;
  *) awk 'NR > 1 && /^set /{exit} NR > 1' "$0"; exit 2 ;;
esac
