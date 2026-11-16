#!/bin/sh
# small_sweep.sh: the phase on small grid instances (roots; nug20 also a depth-3 node), certified bound at marks.
# Rows: the phase's beta factor (wc; 0.36 = today's, wc unset), then projection and extrapolation variants at wc 1.44 / 2.88.
D=/home/user/qap/experiments/bound_lab/cpfast; S=$D/scripts/small.sh; P=$D/par
export BIN=/home/user/wt-cpfast/bin/dnncuts.v4 MARKS="20 50 100 200 300 400 600 800 1200 1600"
for job in "nug12 3:4 - 100" "nug15 3:5 - 100" "nug20 4:5 - 200" "nug20 4:5 3:0,7:6,12:19 200"; do
  set -- $job
  echo "== $1 grid $2 node $3 (first solve $4 iterations, cold)"
  for wc in 0.36 0.72 1.44 2.88 5.76 11.5; do printf "wc %-28s" $wc; SPEC=$3 FIRST=$4 PAR=$P/wc$wc.par $S $1 $2 1600; done
  for wc in 1.44 2.88; do
    for opt in "-cpsweeps 3" "-cpfree 0.05" "-cpfree 0.05 -cpsweeps 2 -cppasses 2" "-cpomega 1.9" "-cpextrap 100:1:30:2" "-cpextrap 100:2:30:2" "-cpextrap 100:2:30:1" "-cpextrap 100:2:30:3"; do
      printf "wc %-4s %-23s" $wc "$opt"; SPEC=$3 FIRST=$4 PAR=$P/wc$wc.par $S $1 $2 1600 $opt; done; done
done
