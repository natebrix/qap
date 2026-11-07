#!/bin/bash
# A second node for the CP-ADMM test: the dual test's child d2 (39:27,30:48,25:49,26:5;
# step-7 Polyak 48,789.87 after the quiet round on the laptop) or any NODE=id:fixes.
# The same two runs as for node 13, dual-test settings:
#   ref   the Polyak protocol (run_ref.sh's flags) with the round-12 pool written
#   cp    CP-ADMM s1p1 from the give-up point with that pool, ITERS (2000) iterations
#   bash run_node.sh ref CORE | cp CORE        logs: polyak_new_<id>.log, cp_s1p1_<id>.log
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../../.." && pwd)
BIN=${BIN:-$R/experiments/dnn_cuts/dnncuts}
NODE=${NODE:-d2:39:27,30:48,25:49,26:5}
id=${NODE%%:*}; x=${NODE#*:}
ITERS=${ITERS:-2000}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$R" || exit 1
case $1 in
  ref) f=$D/logs/polyak_new_$id.log
       grep -q '^FINAL' "$f" 2>/dev/null || taskset -c $2 "$BIN" -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 -i 300 -r 12 -delta 50 -nocap \
         -freeze 0:1500 -vc -x $x -savepoolall 12:$D/logs/pool12_$id.txt qaplib/wil50.dat > "$f.tmp" 2>&1 && mv "$f.tmp" "$f" ;;
  cp)  f=$D/logs/cp_s1p1_$id.log
       grep -q '^FINAL' "$f" 2>/dev/null || taskset -c $2 "$BIN" -grid 5:10 -u 48816 -T 1 -v -I 400 -x $x -cpadmm $D/logs/pool12_$id.txt \
         -cpsweeps 1 -cppasses 1 -cpdiag 200:30 -cpiters $ITERS -v qaplib/wil50.dat > "$f.tmp" 2>&1 && mv "$f.tmp" "$f" ;;
  *) echo "usage: run_node.sh ref|cp CORE"; exit 1 ;;
esac
