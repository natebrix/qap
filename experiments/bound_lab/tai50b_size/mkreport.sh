#!/bin/bash
# Tables for the README from the runs:  ./mkreport.sh   -> logs/*.md, logs/estimate.txt, logs/levels.txt
cd "$(dirname "$0")"
python3 nodes.py runs/root_short/dnn.txt runs/root_long/dnn.txt > logs/roots.md
python3 nodes.py runs/s1_*/dnn.txt > logs/depth1_short.md
python3 nodes.py runs/L1_*/dnn.txt > logs/depth1_long.md
python3 nodes.py runs/spine_*/dnn.txt > logs/spine.md
for s in reh LA0 $(cd runs && ls -d LA_* LB* 2>/dev/null); do
  python3 dives.py runs/$s > logs/dives_$s.txt
  seeds=$(cd runs/$s && for d in s*; do [ -s $d/dives.txt ] && echo ${d#s}; done)
  [ -n "$seeds" ] && python3 -I ../../tai50b_feasibility/estimate/pool.py runs/$s $seeds > logs/pool_$s.txt 2>&1
done
python3 estimate.py > logs/estimate.txt 2>&1; cat logs/estimate.txt
python3 validity.py > logs/validity.txt
python3 webcheck.py > /dev/null; python3 webcost.py > /dev/null
