#!/bin/bash
# depth-4 children of two near-miss depth-3 wil50 nodes: one free facility each, 8 random free locations
cd /Users/nathan_brixius/dev/qap/latest
export VECLIB_MAXIMUM_THREADS=1
OUT=/tmp/claude-501/wil50_d4
python3 - > $OUT/nodes.txt <<'PY'
import random
random.seed(20261030)
parents = {"p9": "26:27,8:38,34:28", "p6": "20:12,40:28,17:25"}
for name, spec in parents.items():
    used_f = {int(a.split(":")[0]) for a in spec.split(",")}
    used_l = {int(a.split(":")[1]) for a in spec.split(",")}
    f = random.choice([x for x in range(50) if x not in used_f])
    locs = random.sample([x for x in range(50) if x not in used_l], 8)
    for l in locs:
        print(f"{name} {spec},{f}:{l}")
PY
cat -n $OUT/nodes.txt | xargs -P 8 -L 1 bash -c 'k=$0; par=$1; spec=$2; ./dnntest -f -T 1 -H 300 -i 800 -u 48816 -x $spec ../qaplib/wil50.dat > /tmp/claude-501/wil50_d4/node$k.log 2>&1'
echo "=== all 16 done ==="
for k in $(seq 1 16); do printf "%-3s %2d %-24s " "$(sed -n ${k}p $OUT/nodes.txt | cut -d' ' -f1)" $k "$(sed -n ${k}p $OUT/nodes.txt | cut -d' ' -f2)"; grep -o "lb=[0-9.]* primal=[0-9.]* status=[a-z_]* iters=[0-9]* cpu=[0-9.]*" $OUT/node$k.log | awk '{split($1,a,"="); printf "lb=%.1f short_by=%.1f %s %s %s\n", a[2], 48814.05-a[2], $3, $4, $5}'; done
