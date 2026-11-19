#!/bin/bash
# tai35b's twin 19:12,18:26 (depth 2) bounded from its parent's state and pool (runs/twin.state, from runs/d1),
# once with the unpatched build and once per mu with the hook (pair 19,18 at 12,26 only); -d 2: bound, branch, stop
H=$(cd "$(dirname "$0")" && pwd); Q=/home/user/wt-cognac/qaplib
one() { $H/run.sh $1 $2 -p $H/tai35b.par -u 283315446 -r $H/runs/twin.state --allow-settings-change -d 2 -o $H/runs/$1/out.state $Q/tai35b.dat; }
[ -f $H/runs/tw_orig/nodes.txt ] || QAP_MIRROR= one tw_orig orig
for mu in ${MUS:-0 0.0625 0.25 1 4 16 0.125 0.5 2 8}; do [ -f $H/runs/tw_mu$mu/nodes.txt ] || QAP_MIRROR="19,18,$mu,12,26" one tw_mu$mu hook; done
