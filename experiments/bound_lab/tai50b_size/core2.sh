#!/bin/bash
# core 2 from 03:25: the twin probe (the path of the best solution found with 31:32, depths 2-4), then B dives again
cd "$(dirname "$0")"
pkill -TERM -f "QAP_EST_SEED=302"; pkill -TERM -f "d2/B.state"; sleep 2
./spine.sh runs/d2/B.state 2 2 4 logs/perm_31_32.txt twin32
./queue.sh 2 "LB 304 runs/d2/B.state 2 6000"
