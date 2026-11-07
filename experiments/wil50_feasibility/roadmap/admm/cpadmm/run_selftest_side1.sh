#!/bin/bash
# The CP self-test on side 1 (README.md, "Side 1", validity (a)), three runs at the
# same time on cores 1-3, logs in logs/side1/:
#   cp_selftest_pool12_13_side1.log  wil50t.dat node 13 (fixes transposed), side 1, the
#                                    round-12 pool of logs/pool12_13.txt with its side field 1
#   cp_selftest_pool12_13_side0.log  wil50.dat node 13, side 0, the same pool
#   cp_selftest_nug12_side1.log      nug12's root, grid 3:4 on its second matrix (side 1),
#                                    its own 3-round pool (rounds_nug12_side1.log)
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../../.." && pwd)
BIN=${BIN:-$R/experiments/dnn_cuts/dnncuts}
L=$D/logs/side1
mkdir -p "$L"
cd "$R" || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
sed 's#^P1/0/#P1/1/#' "$D/logs/pool12_13.txt" > "$L/pool12_13_side1.txt.tmp"
taskset -c 1 "$BIN" -grid 5:10 -side 1 -u 48816 -T 1 -k 1000 -v -I 10 -x 28:21,41:8,19:28 \
  -cpadmm "$L/pool12_13_side1.txt.tmp" -cpiters 10 -cpcheck qaplib/wil50t.dat > "$L/cp_selftest_pool12_13_side1.log" 2>&1 &
taskset -c 2 "$BIN" -grid 5:10 -u 48816 -T 1 -k 1000 -v -I 10 -x 21:28,8:41,28:19 \
  -cpadmm "$D/logs/pool12_13.txt" -cpiters 10 -cpcheck qaplib/wil50.dat > "$L/cp_selftest_pool12_13_side0.log" 2>&1 &
( taskset -c 3 "$BIN" -grid 3:4 -side 1 -c 578 -T 1 -k 150 -v -I 200 -i 200 -r 3 -delta 20 -nocap \
    -savepoolall "3:$L/pool3_nug12_side1.txt" qaplib/nug12.dat > "$L/rounds_nug12_side1.log" 2>&1
  taskset -c 3 "$BIN" -grid 3:4 -side 1 -c 578 -T 1 -v -I 200 -cpadmm "$L/pool3_nug12_side1.txt" -cpiters 800 \
    -cpsweeps 1 -cppasses 1 -cpcheck qaplib/nug12.dat > "$L/cp_selftest_nug12_side1.log" 2>&1 ) &
wait
rm -f "$L/pool12_13_side1.txt.tmp"
grep -H -E "cp-selftest|cp solve" "$L"/cp_selftest_*.log | sed 's#^.*/##'
