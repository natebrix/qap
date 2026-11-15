#!/bin/bash
# after the tai50b nodes: roots of tai40b, tai35b, tai25b, tai20b (plain, comb, top), cutoff 2 x optimum
cd /home/user/qap/experiments/bound_lab/regions
while pgrep -f nodes50.sh > /dev/null; do sleep 5; done
for io in tai40b:1274501896 tai35b:566630890 tai25b:688711292 tai20b:244910638; do
  i=${io%%:*}; c=${io##*:}
  ./run.sh plain_$i $i $c -prefix -I 300 -cpplain -cpiters 2000
  ./run.sh r${i#tai}_comb $i $c -regions regions/$i.comb.txt -side 0 -I 300 -k 1000 -cpsep 100:1000 -cpiters 2000
  ./run.sh r${i#tai}_top $i $c -regions regions/$i.top.txt -side 0 -I 300 -k 1000 -cpsep 100:1000 -cpiters 1000
done
echo TAIL DONE
