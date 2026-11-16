#!/bin/bash
# nodes on the path to the best known solution: the node's optimum is the instance's, so no certified
# value may exceed it. Cutoff = optimum + 1000, so the run goes on and "fathomed" would be a violation.
cd /home/user/qap/experiments/bound_lab/regions
while pgrep -f tail.sh > /dev/null; do sleep 5; done
while read inst val d spec; do
  [ $inst = tai50b ] && continue
  TMO=600 ./run.sh onpath_${inst}_d$d $inst $((val + 1000)) -x $spec -regions regions/$inst.comb.txt -side 0 -I 300 -k 1000 -cpsep 100:1000 -cpiters 2500
done < onpath.txt
echo ONPATH DONE
