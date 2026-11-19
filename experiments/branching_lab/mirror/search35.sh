#!/bin/bash
# the whole tai35b search from the root's state (runs/root/root.state, the root bounded and branched by the unpatched
# build) to the end, one process, with the hook's automatic rule: QAP_MIRROR_AUTO=F:CAP (default 2:64); checkpoints
# every 600 s. ./search35.sh NAME [F:CAP] ; F:CAP = off runs the hook build with no cut (the baseline)
H=$(cd "$(dirname "$0")" && pwd); Q=/home/user/wt-cognac/qaplib; NAME=$1; AUTO=${2:-2:64}
[ "$AUTO" = off ] && AUTO=
QAP_MIRROR_AUTO=$AUTO QAP_MIRROR_CHECK=1 STORE=$H/runs/store_$NAME $H/run.sh $NAME hook3 -p $H/tai35b.par -u 283315446 -r $H/runs/root/root.state --allow-settings-change -k 600 -o $H/runs/$NAME/out.state $Q/tai35b.dat
