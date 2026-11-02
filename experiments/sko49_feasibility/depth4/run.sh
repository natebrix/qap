#!/bin/bash
# depth-4 children (depth4/nodes.txt: "<parent> <spec>") through dnntest, same settings as depth3/run.sh.
# usage: run.sh P [first last]   (lines first..last of nodes.txt; default all)
cd "$(dirname "$0")/../../../latest" || exit 1
OUT=../experiments/sko49_feasibility/depth4
P=${1:-2}; FIRST=${2:-1}; LAST=${3:-1000}
cat -n $OUT/nodes.txt | awk -v a=$FIRST -v b=$LAST '$1>=a && $1<=b' | xargs -P $P -L 1 bash -c 'k=$0; par=$1; spec=$2; stdbuf -oL ./dnntest -T 1 -g 1.0 -G 1.0 -d 4 -H 0 -i 1500 -k 10 -u 23386 -v -x $spec ../qaplib/sko49.dat > '$OUT'/node$k.log 2>&1'
echo "=== done $FIRST-$LAST ==="
