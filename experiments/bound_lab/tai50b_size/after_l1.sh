#!/bin/bash
# After the depth-1 stage (runs/L1_*): the depth-2 open nodes as two state files,
# runs/d2/A.state under the on-path child 31:7 and runs/d2/B.state under the sampled others.
cd "$(dirname "$0")"; mkdir -p runs/d2
python3 splitstate.py -merge runs/d2/all.state runs/L1_0/out.state runs/L1_1/out.state runs/L1_2/out.state runs/L1_3/out.state
python3 splitstate.py runs/d2/all.state runs/d2/A.state -pick 31:7
python3 splitstate.py runs/d2/all.state runs/d2/B.state -drop 31:7
