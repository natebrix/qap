#!/bin/zsh
# usage: run_children.sh f1 f2 ...   (facilities; children = J1 orbit cells of the 5x10 grid)
set -u
S=$(dirname $0)
cd ~/dev/qap/latest
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 DNN_R32_NONCERT=1 DNN_BETA_SCHED=two,f=1.6667,t=60 \
       DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1
echo "== start $(date +%H:%M:%S) facilities $* load $(sysctl -n vm.loadavg)"
for f in "$@"; do for c in 0 1 2 3 4 10 11 12 13 14 20 21 22 23 24; do echo "$f $c"; done; done | \
xargs -P 8 -L 1 sh -c 'f=$0; c=$1; out='"$S"'/children/f${f}_c${c}.log; [ -s "$out" ] && grep -q RESULT "$out" || ./dnntest -T 1 -d 3 -L 2 -r 40 -g 1.0 -G 1.0 -H 0 -i 1500 -k 100 -u 48816 -x $c:$f ../qaplib/wil50.dat > "$out" 2>&1'
echo "== done $(date +%H:%M:%S)"
