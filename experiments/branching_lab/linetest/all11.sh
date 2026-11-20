#!/bin/bash
# all11.sh: rule 2 at every branched depth-17 node of the c47 descent (11; three done earlier in runs/arm17/*_rule2) and
# every branched depth-16 node (7), one level, from the descent's states rebuilt where needed (README). Rebuild chains
# R1-R5 (rebuild.py, tai60b.par = the descent's rule 9), then the arms (armD.sh, rule2.par, latest_dbg), at most 4
# processes; an arm starts when its source state's rebuild has finished.
D=$(cd "$(dirname "$0")"; pwd); cd $D; R2=$D/rule2.par; RB=$D/runs/rebuild
P=44:45,43:58,1:33,6:24,7:19,0:32,10:39
A=$P,46:47,36:53,32:40,3:35; Q=$P,47:47,31:10,46:38,36:53
slot() { while [ $(pgrep -c -x bnbdrv) -ge 4 ]; do sleep 20; done; }
ready() { [ -z "$1" ] || grep -q '^done' $RB/$1.out 2>/dev/null; }
rb() { slot; setsid nohup python3 -u $D/rebuild.py $1 $2 $3 > $RB/$1.out 2>&1 < /dev/null & sleep 15; }
rb R1 $RB/a849/d15.out.state $A,11:36,13:59,12:0,9:55,4:23
rb R2 $RB/a849/d11.out.state $A,11:59,31:10,47:38,51:30,50:6
rb R3 $RB/p2740/d14.out.state $Q,27:40,22:31,11:59,3:35,51:30
rb R4 $RB/path/d13.out.state $Q,32:40,3:35,11:36,13:59,12:0
rb R5 $RB/path/d15.out.state $Q,32:40,3:35,11:59,51:30,50:6
# NAME  wait-for-rebuild  SRC  DEPTH  FIX
cat > $D/runs/all11_jobs.txt <<J
s17_25_31 - $RB/path/d16.out.state 17 25:31
s16_4_23 - $RB/a849/d15.out.state 16 4:23
s16_4_29 - $RB/a849/d15.out.state 16 4:29
s16_48_6 - $RB/path/d15.out.state 16 48:6
s16_50_6 - $RB/path/d15.out.state 16 50:6
s17_8_29 R1 $RB/R1/d16.out.state 17 8:29
s17_9_0 R5 $RB/R5/d16.out.state 17 9:0
s17_9_55q R5 $RB/R5/d16.out.state 17 9:55
s17_9_29 R4 $RB/R4/d16.out.state 17 9:29
s17_9_55p R4 $RB/R4/d16.out.state 17 9:55
s16_12_0 R4 $RB/R4/d15.out.state 16 12:0
s17_48_6 R3 $RB/R3/d16.out.state 17 48:6
s16_51_30c R3 $RB/R3/d15.out.state 16 51:30
s17_9_55a R2 $RB/R2/d16.out.state 17 9:55
s16_50_6a R2 $RB/R2/d15.out.state 16 50:6
J
pending=$(cat $D/runs/all11_jobs.txt)
while [ -n "$pending" ]; do
  left=""
  while read n w s d f; do
    [ -z "$n" ] && continue
    if ready $([ "$w" = - ] || echo $w) && [ $(pgrep -c -x bnbdrv) -lt 4 ]; then
      setsid nohup bash $D/armD.sh $n $s $d $f $R2 > /dev/null 2>&1 < /dev/null & sleep 15
    else left="$left$n $w $s $d $f"$'\n'; fi
  done <<< "$pending"
  pending=$left; [ -n "$pending" ] && sleep 30
done
echo all11 jobs launched
