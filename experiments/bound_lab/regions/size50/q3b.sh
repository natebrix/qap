#!/bin/bash
cd /home/user/qap/experiments/bound_lab/regions/size50
while pgrep -f "q3.sh" > /dev/null; do sleep 5; done
./kid.sh d1_31_23_wc576 kidW576 runs/rootL/root.state 1 31:23
echo Q3B DONE
