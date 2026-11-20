#!/bin/bash
# run1.sh TAG INST FIXES SLACK KICKS[:KICKS3] KICKS2 LINES: t0 (with the kicks3 fix) -> t0out/harvest/TAG.txt
W=/home/user/wt-cognac; H=$(dirname "$0")
/tmp/claude-0/t0h $W/qaplib/$2.dat $W/QAPSOLN/$2.sln "$3" $4 $5 $6 2.5 0 $7 > $H/$1.txt 2>&1
