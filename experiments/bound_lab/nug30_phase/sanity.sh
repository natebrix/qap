#!/bin/bash
# per run: 'New objective value' lines, the final objective, open nodes left in final.state, the solver's node totals
cd "$(dirname "$0")"
for d in runs/*/; do a=$(basename $d); [ -f $d/done ] || continue
  printf "%-8s new-objective lines %s; %s; open nodes %s; %s\n" $a "$(grep -c 'New objective' $d/log.txt)" "$(grep -o 'best objective [0-9]*' $d/log.txt)" "$(grep -c '^node ' $d/final.state)" "$(grep 'Nodes:' $d/log.txt | tail -1 | sed 's/^ *//')"
done
