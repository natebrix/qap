#!/bin/bash
# A long paired run of the AMX kernel's default against one knob setting (default:
# DNN_CORE_AMX_ORDER=1, the interleaved accumulator chains): N alternations A, B, A, B
# of ../../amx/test_amx_gemm.c on CPU 3, so a drift of the host's load hits both.
# Reports min and median of the kernel's own ms per product (Y split + tile loop).
D=$(cd "$(dirname "$0")" && pwd)
CPU=${CPU:-3}
N=${N:-15}
VAR=${VAR:-"DNN_CORE_AMX_ORDER=1"}
K=${K:-2117}; KB=${KBLOCK:-610}
OUT=${OUT:-$D/logs/pair_amx.txt}
echo "# $(date -u +%FT%TZ) K $K kb $KB N $N variant [$VAR] load $(cut -d' ' -f1-3 /proc/loadavg)" >> "$OUT"
for i in $(seq 1 $N); do
  a=$(env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 taskset -c $CPU "$D/test_amx" $K $KB | grep -o 'gemm [0-9.]* ms' | awk '{print $2}')
  b=$(env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 $VAR taskset -c $CPU "$D/test_amx" $K $KB | grep -o 'gemm [0-9.]* ms' | awk '{print $2}')
  echo "pair $i base $a variant $b" >> "$OUT"
done
python3 - "$OUT" <<'EOF'
import sys, statistics
a=[]; b=[]
for line in open(sys.argv[1]):
    if line.startswith('pair'):
        f=line.split(); a.append(float(f[3])); b.append(float(f[5]))
# the last block only: pairs after the last header
blk=[]
for line in open(sys.argv[1]):
    if line.startswith('#'): blk=[]
    elif line.startswith('pair'): blk.append(line.split())
a=[float(f[3]) for f in blk]; b=[float(f[5]) for f in blk]
d=[y-x for x,y in zip(a,b)]
print("pairs %d: base min %.1f median %.1f | variant min %.1f median %.1f | paired difference median %+.1f ms, variant faster in %d of %d" % (
  len(a), min(a), statistics.median(a), min(b), statistics.median(b), statistics.median(d), sum(1 for x in d if x<0), len(d)))
EOF
