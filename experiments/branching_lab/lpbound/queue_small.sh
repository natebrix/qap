#!/bin/bash
cd "$(dirname "$0")"
for a in "nug12 -" "tai12b -" "had14 0:7" "nug15 0:4" "tai15b 0:0" "tai15b 0:3" "nug15 0:4,1:2"; do
  python3 tsmall.py $a 4000 2>&1 | grep -v "^Running\|^Includes"
done
