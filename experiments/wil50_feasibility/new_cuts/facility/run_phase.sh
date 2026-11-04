#!/bin/bash
# Triangle phase (step-7 protocol: 400 uncut, ROUNDS (8; reference 12) rounds of 300, 1,000 regional
# triangles a round, Polyak toward best + 50 uncapped), stopped before the quiet
# round; saves runs/<id>/tri.* (term, state, round-0 iterate). Then the quiet
# round (1,500 iterations from the best round's multipliers) as a separate
# process, dumping the quiet iterate runs/<id>/quiet.Y.bin.
#   bash run_phase.sh <id> <spec>      (resumable: FINAL lines)
set -u
D=$(cd "$(dirname "$0")" && pwd)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_AMX=0 DNN_CHOLQR_PASSES=1
id=$1; spec=$2
F="$D/src/fdrv -x $spec -u 48816 -d 3 -L 2 -r 40"
mkdir -p $D/runs/$id; cd $D/runs/$id
if ! grep -q '^FINAL' phase.log 2>/dev/null; then
  taskset -c 3 $F -grid 5:10 -k 1000 -I 400 -i 300 -rounds ${ROUNDS:-8} -delta 50 -nocap -freeze 1500 -phase tri /home/user/qap/qaplib/wil50.dat > phase.log.tmp 2>&1 && mv phase.log.tmp phase.log
fi
grep -q '^FINAL' phase.log || exit 1
if ! grep -q '^FINAL' quiet.log 2>/dev/null; then
  taskset -c 3 $F -cont tri -iters 1500 -dumpY quiet -save quiet /home/user/qap/qaplib/wil50.dat > quiet.log.tmp 2>&1 && mv quiet.log.tmp quiet.log
fi
grep '^PHASE\|^RESULT' phase.log quiet.log
