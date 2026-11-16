#!/bin/bash
# q1c.sh: after the chain to depth 8: on down the path to depth 12
D=$(cd $(dirname $0); pwd); cd $D
while pgrep -f "q1b.sh" > /dev/null; do sleep 20; done
prev=$(ls -d runs/d8_*path | head -1 | xargs basename)
./chain.sh $prev 9 12
