#!/bin/bash
# cached.sh MODE [NODE]: how many child solves sb.early = MODE needs at wil50's root or a depth-1 node, at no
# cost in child solves: the firing of ../wil50_sb solved every child (runs/*.progress.txt there), and its
# progress file is reused here under this run's settings hash (the hash covers the params file, which differs
# by its sb.early line only). The node's own solve and cut phase run again (about 20 minutes at depth 1).
set -uo pipefail
D=$(cd "$(dirname "$0")"; pwd); mode=$1; node=${2:-}; name=${node:-root}; name=${name//:/-}; tag=$name.e$mode
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$D/../../latest"
args=(-p "$D/par/sb_e$mode.par" -s SYM/wil50t.sym -u 48818 ${node:+-F "$node"} -d $([ -n "$node" ] && echo 1 || echo 0) ../qaplib/wil50t.dat)
new=$(timeout 20 ./bnbdrv "${args[@]}" -o "$D/runs/$tag.tmp.state" 2>&1 | awk '/^settings [0-9]+ hash/ {print $4; exit}')
old=$(awk 'NR == 1 {print $2}' "$D/../wil50_sb/runs/$name.progress.txt")
[ -n "$new" ] && [ -n "$old" ] || { echo "no hash: new '$new' old '$old'"; exit 1; }
sed "s/^sb-child $old /sb-child $new /" "$D/../wil50_sb/runs/$name.progress.txt" > "$D/runs/$tag.progress.txt"
rm -f "$D/runs/$tag.sb.txt" "$D/runs/$tag.tmp.state"
QAP_SB_TRACE="$D/runs/$tag.sb.txt" QAP_SB_PROGRESS="$D/runs/$tag.progress.txt" ./bnbdrv "${args[@]}" -o "$D/runs/$tag.state" > "$D/runs/$tag.out.txt" 2>&1
echo "$tag exit $? $(date -u +%T)" >> "$D/runs/done.txt"
