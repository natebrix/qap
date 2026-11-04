#!/bin/bash
# after queue f: scale of the DNN's 3-local inconsistency, nug12 / nug20 / nug30 roots
V=/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/venv/bin/python
cd /home/user/qap/experiments/wil50_feasibility/new_cuts/literature/small
until grep -q "QUEUE DONE" logs/run_queue_f.out 2>/dev/null; do sleep 20; done
T=/home/user/qap/latest/dnntest
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for inst in nug12 nug20 nug30; do
  [ -s logs/${inst}_root.Y.bin ] || taskset -c 1 $T -u 1e9 -i 3000 -H 0 -dumpY logs/${inst}_root /home/user/qap/qaplib/$inst.dat > logs/${inst}_root.log 2>&1
  grep -q RESULT logs/${inst}_lmpscale.log 2>/dev/null || taskset -c 1 $V -I lmp_scale.py /home/user/qap/qaplib/$inst.dat logs/${inst}_root --heavy 30 --rand 30 > logs/${inst}_lmpscale.log 2>&1
done
echo QUEUE G DONE
