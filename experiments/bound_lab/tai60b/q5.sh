#!/bin/bash
# q5.sh: second slot after q4: phases of the two children of the survivor 1:39 that its first-solve sample left
# open (22:48, 0.62M short; 22:9, 0.24M short; both facilities of location 22's cluster in the best known solution)
D=$(cd $(dirname $0); pwd); cd $D
while pgrep -f "q4.sh" > /dev/null; do sleep 10; done
./kid.sh d5_1_39_22_48 kidW runs/d4_1_39/out.state 5 22:48
./kid.sh d5_1_39_22_9 kidW runs/d4_1_39/out.state 5 22:9
