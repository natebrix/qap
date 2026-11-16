#!/bin/bash
# One whole nug30 search, one arm, 4 workers (bench/parallel.py):  ./run.sh ARM [JOBS]
#   -> runs/ARM/{log.txt,nodes.txt,dnn.txt,time.txt,cmd.txt,final.state}
set -u
H=$(cd "$(dirname "$0")" && pwd); R=$(cd "$H/../../.." && pwd)
ARM=$1; JOBS=${2:-4}; UB=${UB:-6124}; PAR=${PAR:-$ARM}; D=$H/runs/$ARM; mkdir -p "$D"
[ -f "$D/done" ] && { echo "$D done"; exit 0; }
rm -f "$D/nodes.txt" "$D/dnn.txt"; W=$(mktemp -d ${TMPDIR:-/tmp}/n30-$ARM-XXXX)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export QAP_NODE_TRACE=$D/nodes.txt QAP_DNN_TRACE=$D/dnn.txt QAP_DNN_TRACE_HIST=1
cd "$R"
CMD="python3 bench/parallel.py qaplib/nug30.dat --jobs $JOBS --ub $UB --params $H/par/$PAR.par --workdir $W --checkpoint $D/final.state --checkpoint-every 100000"
{ echo "$CMD"; git rev-parse --short HEAD; date -u +%FT%TZ; uptime; } > "$D/cmd.txt"
TIMEFORMAT="wall %R user %U sys %S"; { time $CMD > "$D/log.txt" 2> "$D/err.txt"; } 2> "$D/time.txt"
{ echo "exit $?"; date -u +%FT%TZ; uptime; } >> "$D/cmd.txt"
rm -rf "$W"; touch "$D/done"
