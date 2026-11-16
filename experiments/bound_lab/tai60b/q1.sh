#!/bin/bash
# q1.sh: first tai60b slot: the on-path chain to depth 4 (depth 2 was started by hand), then first solves at
# depths 3 and 4 under the on-path nodes
D=$(cd $(dirname $0); pwd); cd $D
until grep -q "^exit" runs/d2_43_58path/cmd.txt; do sleep 20; done
./chain.sh d2_43_58path 3 4
