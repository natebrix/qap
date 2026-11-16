#!/bin/bash
# One complete subtree, one arm, one solver process (cold from -F):
#   ./run.sh NAME ARM SPEC [EXTRA]     -> runs/NAME_ARM/{log.txt,nodes.txt,dnn.txt,cmd.txt}
# EXTRA: more bnbdrv flags, e.g. '-d 3' (bound the root only: a probe).
# ARM names par/ARM.par. Instance, symmetry file and incumbent as experiments/estimates/tho40/settings.sh.
# NAME path*: incumbent 240,518 (the path check). NAME w* (not wroot): a WARM probe of the node SPEC: its parent's
# state runs/wroot_ARM/root.state (written by the job "wroot ARM 26:10,8:0 -d 2 -o .../root.state"), cut to this
# child (cutstate.py), is resumed with -r, so the node starts from its parent's shipped DNN state as in a run.
# (Replace this file only by mv of a new one: bash reads a running script from the file.)
set -u
H=$(cd "$(dirname "$0")/.." && pwd); R=/home/user/qap
# stop/run2.sh: run.sh for the stop-rule study. Arms named *h* (the cut.handover key) run the patched solver of
# the worktree /home/user/wt-tho40search (stop/handover.patch); every other arm runs latest/bnbdrv of the main checkout.
case $2 in *h*) R=/home/user/wt-tho40search;; esac
# ARM may be written ARM@ROOTARM: a warm probe then resumes ROOTARM's root state under ARM's settings
# (--allow-settings-change; used where the two arms differ in nothing the root's bound depends on).
ROOTARM=${2##*@}; set -- "$1" "${2%%@*}" "${@:3}"
NAME=$1; ARM=$2; SPEC=$3; EXTRA=${4:-}; UB=${UB:-240516}
case "$SPEC" in *" "*) EXTRA="${SPEC#* } $EXTRA"; SPEC="${SPEC%% *}";; esac   # a caller that passed "SPEC EXTRA" as one word
case $NAME in path*) UB=240518;; esac   # the path check: the best known + 2, as tho40-b and path_check/dive.py
D=$H/runs/${NAME}_$ARM; mkdir -p "$D"
[ -f "$D/done" ] && { echo "$D already done"; exit 0; }
rm -f "$D/nodes.txt" "$D/dnn.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt QAP_DNN_TRACE_HIST=1
START="-F $SPEC"
case $NAME in
  wroot*) ;;
  w*) NF=$(echo "$SPEC" | tr ',' '\n' | wc -l)
      python3 "$H/cutstate.py" "$H/runs/wroot_$ROOTARM/root.state" "$D/in.state" "$NF" "${SPEC##*,}" || { echo "cutstate failed for $NAME"; exit 1; }
      START="-r $D/in.state"; [ "$ROOTARM" = "$ARM" ] || START="$START --allow-settings-change";;
esac
cd $R/latest
CMD="./bnbdrv -p $H/par/$ARM.par -s SYM/tho40t.sym -u $UB $START $EXTRA /home/user/qap/qaplib/tho40t.dat"
{ echo "$CMD"; echo "QAP_NODE_TRACE=$QAP_NODE_TRACE QAP_DNN_TRACE=$QAP_DNN_TRACE QAP_DNN_TRACE_HIST=1 OPENBLAS_NUM_THREADS=1"; date -u +%FT%TZ; uptime; } > "$D/cmd.txt"
$CMD > "$D/log.txt" 2> "$D/err.txt"
echo "exit $?" >> "$D/cmd.txt"; date -u +%FT%TZ >> "$D/cmd.txt"; uptime >> "$D/cmd.txt"
touch "$D/done"
