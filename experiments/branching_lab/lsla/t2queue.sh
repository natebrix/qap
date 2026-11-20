#!/bin/bash
# the tai35b arms of T2, one after another (2 slots each)
D=$(cd "$(dirname "$0")"; pwd)
P=${1:-tai35b}
# an arm whose cmd.txt ends "exit 0" has finished and is skipped (a restart reruns the rest)
arm() { grep -q "^exit 0 " "$D/t2/${P}_$1/cmd.txt" 2>/dev/null && return 0; rm -rf "$D/t2/${P}_$1"; "$D/t2.sh" "$P" "$@"; }
arm prod
arm a QAP_LS=a:2:0:2
arm b2 QAP_LS=b:2:0:2 QAP_LS_CAP=30
arm b4 QAP_LS=b:4:0:2 QAP_LS_CAP=30
