#!/bin/bash
# search.sh NAME PARFILE SPEC [EXTRA]: as ../tho40_search/run.sh with the worktree's bnbdrv.
#   NAME wroot*: cold from -F SPEC (use EXTRA '-d 2 -o search/NAME/root.state');
#   NAME w*: the warm probe of node SPEC from search/wroot_<par name>/root.state (cutstate.py), resumed with -r.
set -u
D0=/home/user/qap/experiments/bound_lab/cpfast; W=/home/user/wt-cpfast; H=/home/user/qap/experiments/bound_lab/tho40_search
NAME=$1; PAR=$2; SPEC=$3; EXTRA=${4:-}; UB=${UB:-240516}
ARM=$(basename $PAR .par)
D=$D0/search/${NAME}; mkdir -p "$D"
rm -f "$D/nodes.txt" "$D/dnn.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt QAP_DNN_TRACE_HIST=1
START="-F $SPEC"
case $NAME in
  wroot*) ;;
  w*) NF=$(echo "$SPEC" | tr ',' '\n' | wc -l)
      python3 "$H/cutstate.py" "${WROOT:-$D0/search/wroot_$ARM/root.state}" "$D/in.state" "$NF" "${SPEC##*,}" || { echo "cutstate failed for $NAME"; exit 1; }
      START="-r $D/in.state";;
esac
cd $W/latest
CMD="${BNB:-$W/bin/bnbdrv.v3} -p $PAR -s SYM/tho40t.sym -u $UB $START $EXTRA ../qaplib/tho40t.dat"
{ echo "$CMD"; date -u +%FT%TZ; uptime; } > "$D/cmd.txt"
$CMD > "$D/log.txt" 2> "$D/err.txt"
echo "exit $?" >> "$D/cmd.txt"; date -u +%FT%TZ >> "$D/cmd.txt"; uptime >> "$D/cmd.txt"
rm -f "$D/in.state"
touch "$D/done"
