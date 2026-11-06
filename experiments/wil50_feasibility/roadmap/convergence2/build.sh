#!/bin/bash
# Build dnntest, dnncuts and bnbdrv in this worktree at low priority on core 0
# (the box is shared; see REPORT.md, compute rules). A checkout from before
# Nov 6 still has src/MESCHACH, which its solver links: that is built first.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../.." && pwd)
if [ -f "$R/src/MESCHACH/makefile" ]; then
  cd "$R/src/MESCHACH" && nice -n 19 taskset -c 0 make basic > "$D/logs/build_meschach.log" 2>&1 || exit 1
fi
cd "$R" && nice -n 19 taskset -c 0 make -C latest dnntest dnncuts bnb > "$D/logs/build_latest.log" 2>&1
echo "build exit $?"
