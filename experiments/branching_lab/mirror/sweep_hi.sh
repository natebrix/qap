#!/bin/bash
# the same twin, cutoff raised to 290,000,000 (incumbent line of the state file edited) so that the solves do not stop at the real
# cutoff: the bound of cost + mu Delta beyond the mu that fathoms, to find the peak of the concave curve
H=$(cd "$(dirname "$0")" && pwd); Q=/home/user/wt-cognac/qaplib
one() { $H/run.sh $1 $2 -p $H/tai35b.par -u 290000000 -r $H/runs/twin_hi.state --allow-settings-change -d 2 -o $H/runs/$1/out.state $Q/tai35b.dat; }
for mu in ${MUS:-0 8 16 32 64}; do [ -f $H/runs/hi_mu$mu/nodes.txt ] || QAP_MIRROR="19,18,$mu,12,26" one hi_mu$mu hook; done
