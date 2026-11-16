#!/bin/bash
# the rest of the queue (replaces the end of tail.sh and onpath.sh): tai35b/25b/20b roots, on-path nodes,
# a tai50b on-path node (15 minutes), then a search of tai30b with the cut phase (20 minutes)
cd /home/user/qap/experiments/bound_lab/regions
until grep -q "r40b_comb.log:FINAL" logs/tail.out; do sleep 3; done
pkill -f tail.sh; sleep 0.5; pkill -f "wt-regions/experiments/dnn_cuts/dnncuts"; sleep 1
for io in tai35b:566630890 tai25b:688711292 tai20b:244910638; do
  i=${io%%:*}; c=${io##*:}
  ./run.sh plain_$i $i $c -prefix -I 300 -cpplain -cpiters 2000
  ./run.sh r${i#tai}_comb $i $c -regions regions/$i.comb.txt -side 0 -I 300 -k 1000 -cpsep 100:1000 -cpiters 2000
done
echo ROOTS DONE
# nodes on the path to the best known solution: the node's optimum is the instance's, so no certified value
# may exceed it. Cutoff = optimum + 1000: the run goes on, and "fathomed" would be a violation.
while read inst val d spec; do
  [ $inst = tai30b ] || continue
  TMO=300 ./run.sh onpath_${inst}_d$d $inst $((val + 1000)) -x $spec -regions regions/$inst.comb.txt -side 0 -I 300 -k 1000 -cpsep 100:1000 -cpiters 2500
done < onpath.txt
TMO=900 ./run.sh onpath_tai50b_d3 tai50b 458822517 -x 7:40,40:39,42:43 -regions regions/tai50b.comb.txt -side 0 -I 300 -k 1000 -cpsep 100:1000 -cpiters 3000
echo ONPATH DONE
cd /home/user/wt-regions/latest
timeout 1200 ./bnbdrv -p /home/user/qap/experiments/bound_lab/regions/bnb/tai30b_cut.par -u 637117114 ../qaplib/tai30b.dat > /home/user/qap/experiments/bound_lab/regions/logs/bnb_tai30b_cut.log 2>&1
echo BNB DONE $?
