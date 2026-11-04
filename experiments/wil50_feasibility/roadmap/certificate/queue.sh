#!/bin/bash
# The measurement queue of REPORT.md, one line of work per core (cores 2 and
# 3 only, one process each, single-threaded). Dumps (39 MB each) go to $DUMP,
# outside the repository.
#   queue.sh 2   E1: cold children of node 10 that fathom slowly in the search:
#                47:45 (450 its, also with the line search), the depth-5
#                49:45,0:12 (650 its), 37:45, 27:45, 7:45, 6:45, 46:45; then
#                E1w: 7:45 and 47:45 warm-started from node 10's state
#   queue.sh 3   E4 (node 10 dump, child bounds on column 45), E2 on [9:45]
#                (phase miss, dump, child bounds on row 19), E2 on [48:45]
D=$(cd "$(dirname "$0")" && pwd)
DUMP=${DUMP:-/tmp/certx_dumps}
mkdir -p "$DUMP"
case $1 in
  2)
    "$D/run_e1.sh" 2 175 47
    "$D/run_e1.sh" 2 171 49:45,0:12 7 37
    "$D/run_e1w.sh" 2 171 7 47
    "$D/run_e1.sh" 2 171 27 6 46
    ;;
  3)
    "$D/run_e2.sh" 3 171 n48_45 "39:27,30:48,25:49,48:45" 470 10 12 "$DUMP"
    "$D/run_e4.sh" 3 "$DUMP"
    "$D/run_e2.sh" 3 187 n9_45 "39:27,30:48,25:49,9:45" 430 10 12 "$DUMP"
    taskset -c 3 python3 "$D/scripts/child_gl.py" "$DUMP/n9_45" line:row:18 > "$D/logs/e2_n9_45_childgl_row19.log" 2>&1
    ;;
esac
