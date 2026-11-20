#!/bin/bash
# scan.sh: the rule 2 / rule 9 scan at depths 11, 13, 15 (README). The three branched nodes per depth with the most
# branching children under rule 9's line (out/scan_select.txt), rebuilt where needed (rebuild.py) and bounded under
# rule 2 one level (armD.sh, rule2.par). JOBS processes at a time (default 4).
D=$(cd "$(dirname "$0")"; pwd); cd $D; R2=$D/rule2.par; RB=$D/runs/rebuild
P=44:45,43:58,1:33,6:24,7:19,0:32,10:39
JOBS=${JOBS:-4}
wait_slot() { while [ $(pgrep -c -x bnbdrv) -ge $JOBS ]; do sleep 20; done; }
job() { wait_slot; setsid nohup bash -c "$*" > /dev/null 2>&1 < /dev/null & sleep 15; }
cd $RB
job "bash $D/armD.sh d11_a $RB/a849/d10.out.state 11 3:35 $R2"
job "bash $D/armD.sh d11_c $RB/path/d10.out.state 11 36:53 $R2"
job "bash $D/armD.sh d13_a $RB/path/d12.out.state 13 3:35 $R2"
job "bash $D/armD.sh d13_c $RB/a849/d12.out.state 13 13:59 $R2"
job "python3 -u $D/rebuild.py b3440 $D/runs/d9_46_47/out.state $P,46:47,36:53,34:40 > $RB/b3440.out 2>&1 && bash $D/armD.sh d11_b $RB/b3440/d10.out.state 11 3:35 $R2"
job "bash $D/armD.sh d15_a $RB/a849/d14.out.state 15 9:55 $R2"
job "bash $D/armD.sh d15_b $RB/path/d14.out.state 15 51:30 $R2"
job "python3 -u $D/rebuild.py p3440 $RB/path/d11.out.state $P,47:47,31:10,46:38,36:53,34:40 > $RB/p3440.out 2>&1 && bash $D/armD.sh d13_b $RB/p3440/d12.out.state 13 3:35 $R2"
job "python3 -u $D/rebuild.py p2740 $RB/path/d11.out.state $P,47:47,31:10,46:38,36:53,27:40,22:31,11:59 > $RB/p2740.out 2>&1 && bash $D/armD.sh d15_c $RB/p2740/d14.out.state 15 3:35 $R2"
