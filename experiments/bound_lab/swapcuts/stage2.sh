#!/bin/sh
# stage 2: exact ceiling, one process at a time. ./stage2.sh
cd "$(dirname "$0")"
for inst in nug12 had12 scr12 chr12a tai12a tai12b rou12 nug14; do
  ./exch ../../../qaplib/$inst.dat 1 4 5 0.05 > logs/stage2_$inst.out 2>&1
done
