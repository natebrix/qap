#!/bin/bash
# after the long on-path run: the one depth-5 child of the survivor 34:28 that did not fathom at its first solve
# (42:35, 475k short): does a phase (wc = 2.88) fathom it? Needs the survivor's out.state.
cd /home/user/qap/experiments/bound_lab/regions/size50
while pgrep -f "q10.sh" > /dev/null; do sleep 3; done
./kid.sh d5_42_35 kidW288 runs/d4_34_28/out.state 5 42:35
echo Q11 DONE
