#!/bin/bash
# Count-split disjunction on a node with a saved quiet state: pick the split on
# the quiet iterate, write the side terms, append the cont jobs to jobs.txt.
#   split_jobs.sh <id> <spec>
set -e
D=$(cd "$(dirname "$0")/.." && pwd); id=$1; spec=$2
cd $D/scripts
taskset -c 3 python3 pick_split.py $spec $D/runs/$id/quiet.Y.bin > $D/runs/$id/split_pick.log
F="$D/src/fdrv -u 48816 -d 3 -L 2 -r 40 -x $spec -cont quiet -iters 400"
echo "runs/$id split_ctl.log $F /home/user/qap/qaplib/wil50.dat" >> $D/jobs.txt
for side in le ge; do for mu in 10 40 160; do
  ex=$(python3 count_term.py $spec --pick $D/runs/$id/split_pick.log --side $side --mu $mu --out $D/runs/$id/split_${side}_$mu.bin)
  echo "runs/$id split_${side}_$mu.log $F -extra $ex /home/user/qap/qaplib/wil50.dat" >> $D/jobs.txt
done; done
cat $D/runs/$id/split_pick.log
