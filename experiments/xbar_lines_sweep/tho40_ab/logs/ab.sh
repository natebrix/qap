#!/bin/bash
cd ~/dev/qap/latest
export VECLIB_MAXIMUM_THREADS=1
S=/private/tmp/claude-501/-Users-nathan-brixius-dev-brixius-web/f833c398-2e99-4f2e-a75b-8111d7ca8792/scratchpad/tho40_ab
echo "$(date +%H:%M:%S) start"
for node in 26:3,31:7,37:0 25:32,26:8,39:0 24:7,26:1,39:0; do
  for k in 1 2; do
    tag=$(echo $node | tr ':,' '-_')_x$k
    ( /usr/bin/time -p ./bnbdrv -p $S/x$k.par -u 240516 -F $node ../qaplib/tho40t.dat > $S/$tag.log 2> $S/$tag.time; echo "$(date +%H:%M:%S) done $tag: $(grep -E 'Nodes:' $S/$tag.log | tr -s ' ') $(grep -E 'Total Time' $S/$tag.log | tr -s ' ') $(grep real $S/$tag.time)" ) &
  done
done
wait
echo "$(date +%H:%M:%S) all done"
