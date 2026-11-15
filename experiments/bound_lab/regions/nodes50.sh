#!/bin/bash
# tai50b depth-3 nodes of tai50b_feasibility/NOTES.md section 2, true cutoff, 25-minute cap each
cd /home/user/qap/experiments/bound_lab/regions
for nd in "7 0:44,49:49,25:31" "1 20:19,28:24,14:35" "10 5:40,40:22,44:3"; do
  set -- $nd
  TMO=1500 ./run.sh t50_d3_node$1_comb tai50b 458821516.05 -x $2 -regions regions/tai50b.comb.txt -side 0 -I 300 -k 1000 -cpsep 100:1000 -cpiters 3000
done
echo NODES50 DONE
