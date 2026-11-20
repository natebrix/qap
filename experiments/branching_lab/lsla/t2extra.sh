#!/bin/bash
# extra tai35b arms, after t2queue.sh: mode (b) K=2 with the signal's 30 kicks only (no refinement), which at
# 18:26 19:12 ranks c23 (c's choice) first where the refined signal leaves it out of the top 12
D=$(cd "$(dirname "$0")"; pwd)
P=${1:-tai35b}
while pgrep -f "t2queue.sh $P" > /dev/null; do sleep 60; done
arm() { grep -q "^exit 0 " "$D/t2/${P}_$1/cmd.txt" 2>/dev/null && return 0; rm -rf "$D/t2/${P}_$1"; "$D/t2.sh" "$P" "$@"; }
SEARCH=30:0:2.5:300 arm b2r0 QAP_LS=b:2:0:2 QAP_LS_CAP=30
