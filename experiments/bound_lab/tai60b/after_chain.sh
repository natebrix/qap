#!/bin/bash
# after_chain.sh: when the chain below depth 18 ends, build the strata of its levels, queue their dives, and
# start a fourth worker on the freed core
H=$(cd "$(dirname "$0")" && pwd); cd $H
while pgrep -f "chain.sh d18_19_13path" > /dev/null; do sleep 30; done
export OPENBLAS_NUM_THREADS=1
for D in $(seq 20 43); do prev=$(ls -d runs/d$((D-1))_*path 2>/dev/null | grep -v under | head -1); [ -n "$prev" ] && [ -s $prev/out.state ] || break
  python3 -I strata.py $(basename $prev) $D >> logs/strata20.out; last=$D; done
tmp=$(mktemp); for D in $(seq 20 ${last:-19}); do echo "far$D $((100+D)) off${D}_far.state $D 3"; grep -q "^nodes 0" runs/strata/off${D}_near.state || echo "near$D $((200+D)) off${D}_near.state $D 3"; done > $tmp
( flock 9; cat $tmp jobs.txt > jobs.new; mv jobs.new jobs.txt ) 9> runs/jobs.lock
./worker.sh 4
