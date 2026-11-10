#!/bin/bash
cd ~/dev/qap
export VECLIB_MAXIMUM_THREADS=1
S=/private/tmp/claude-501/-Users-nathan-brixius-dev-brixius-web/f833c398-2e99-4f2e-a75b-8111d7ca8792/scratchpad/xbar_sweep
run() { # instance ub par label
  echo "$(date +%H:%M:%S) start $4"
  python3 bench/parallel.py qaplib/$1.dat --jobs 8 --ub $2 --params $S/$3 > $S/$4.log 2>&1
  echo "$(date +%H:%M:%S) done  $4: $(grep -E 'Nodes:' $S/$4.log | tr -s ' ') $(grep -E 'fw iterations' $S/$4.log | tr -s ' ') $(grep -E 'Total Time' $S/$4.log | tr -s ' ') $(grep -oE 'wall [0-9.]+ s' $S/$4.log)"
}
for k in 1 2 4 8; do run tho30 149936 tho30_r9_x$k.par tho30_r9_x$k; done
for k in 1 4 8; do run nug30 6124 nug30_x$k.par nug30_x$k; done
echo "$(date +%H:%M:%S) sweep finished"
