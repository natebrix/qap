#!/bin/bash
# cut phase (cuts/run_cuts.sh protocol) on node 2's depth-4 children, lines 9-14 of depth4/nodes.txt, 4 at a time
cd "$(dirname "$0")"
awk 'NR>=9 && NR<=14 {print NR, $2}' ../depth4/nodes.txt | xargs -P 4 -L 1 bash -c './run_cuts.sh d4_$0 $1'
echo "=== done $(date)"
