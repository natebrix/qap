#!/bin/bash
# regional cuts on a depth-3 miss: the wil50 step-6/7 protocol from the give-up point (400 uncut
# iterations, 12 rounds of 300 with 1,000 regional triangles each, uncapped Polyak toward best + 50,
# quiet round of 1,500 from the best round's multipliers). usage: run_cuts.sh <tag> <spec> [-five 500]
cd "$(dirname "$0")/../../../latest" || exit 1
OUT=../experiments/sko49_feasibility/cuts
tag=$1; spec=$2; shift 2
stdbuf -oL ../experiments/dnn_cuts/dnncuts -x $spec -u 23386 -grid 7:7 -side 0 -T 1 -k 1000 -v -I 400 -i 300 -r 12 -delta 50 -nocap -freeze 0:1500 "$@" ../qaplib/sko49.dat > $OUT/$tag.log 2>&1
