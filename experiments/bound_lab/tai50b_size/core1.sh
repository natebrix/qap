#!/bin/bash
# core 1 from 03:36: the odd spine levels' dive sets, then the tail (waits for the states core0.sh builds)
cd "$(dirname "$0")"
./queue.sh 1 "LA_off5 215 runs/d2/A_off5.state 5 900"
until [ -s runs/d2/A_off7.state ]; do sleep 5; done
./queue.sh 1 "LA_off7 217 runs/d2/A_off7.state 7 900"
until [ -s runs/d2/A_tail8.state ]; do sleep 5; done
./queue.sh 1 "LA_tail8 229 runs/d2/A_tail8.state 9 6000"
