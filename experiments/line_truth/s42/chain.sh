#!/bin/bash
# s42/chain.sh: once d2/kids.sh has finished, the firings at sko42-261124's eight heaviest depth-2 nodes, four at a time.
D=$(cd "$(dirname "$0")"; pwd)
while ! grep -q "^all exit" "$D/../d2/kids/done.txt" 2>/dev/null; do sleep 60; done
printf '%s\n' 22:0,24:41 22:0,10:41 22:1,21:0 22:7,21:0 22:0,28:41 22:7,15:0 22:1,15:0 22:0,40:41 | xargs -P 4 -L 1 "$D/run.sh"
echo "all exit 0 $(date -u +%FT%TZ)" >> "$D/runs/done.txt"
