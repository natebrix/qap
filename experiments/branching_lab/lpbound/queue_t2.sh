#!/bin/bash
# the remaining tai60b truth nodes, PDLP 1e-4, 3 threads
cd "$(dirname "$0")"
P7=44:45,43:58,1:33,6:24,7:19,0:32,10:39
for par in $P7,46:47,36:53,32:40,3:35,11:36,13:59,12:0,9:55,4:29,8:49 $P7,46:47,36:53,32:40,3:35,11:36,13:59,12:0,9:55,4:29,8:57 $P7,14:54 $P7; do
  THREADS=3 python3 t2.py $par 1e-4 2>&1 | grep -v "^Running\|^Includes"
done
P17=44:45,43:58,1:33,6:24,7:19,0:32,10:39,47:47,31:10,46:38,36:53,32:40,3:35,11:59,51:30,48:6,22:31
THREADS=3 python3 childlp.py $P17 c26 6 2>&1 | grep -v "^Running\|^Includes"
