#!/bin/bash
# family sweep at the tai30b root: 300 uncut iterations, then 1,000 CP iterations
cd /home/user/qap/experiments/bound_lab/regions
for f in line2:0 hier:0 top:0 ball4:0 random30:0 linegap40.03:0 fhier:1 fhub10:1 fspec2:1 line4:0; do
  fam=${f%%:*}; side=${f##*:}
  ./run.sh r30_$fam tai30b 1274234226 -regions regions/tai30b.$fam.txt -side $side -I 300 -k 1000 -cpsep 100:1000 -cpiters 1000
done
echo SWEEP30 DONE
