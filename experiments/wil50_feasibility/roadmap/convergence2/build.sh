#!/bin/bash
# Build meschach, dnntest, dnncuts and bnbdrv in this worktree at low priority
# on core 0 (the box is shared; see REPORT.md, compute rules).
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
cd "$R/src/MESCHACH" && nice -n 19 taskset -c 0 make basic > "$D/logs/build_meschach.log" 2>&1 || exit 1
cd "$R" && nice -n 19 taskset -c 0 make -C latest dnntest dnncuts bnb > "$D/logs/build_latest.log" 2>&1
echo "build exit $?"
