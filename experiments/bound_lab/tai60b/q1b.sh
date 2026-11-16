#!/bin/bash
# q1b.sh: after the chain to depth 4: go on down the path to depth 8
D=$(cd $(dirname $0); pwd); cd $D
while pgrep -f "chain.sh d2_43_58path" > /dev/null; do sleep 20; done
prev=$(ls -d runs/d4_*path | head -1 | xargs basename)
./chain.sh $prev 5 8
