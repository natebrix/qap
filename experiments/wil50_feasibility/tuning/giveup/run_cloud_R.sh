#!/bin/bash
# Gap-ratio give-up rule (dnntest -R) on the 32 nodes of all32.txt: -H 150
# -R 0.65, otherwise as run_cloud.sh. Resumable. Logs: cloud/H150R065_<id>.log
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
one() {
  k=$1; spec=$2
  f=$D/cloud/H150R065_$k.log
  grep -q '^RESULT' "$f" 2>/dev/null && exit 0
  "$R/latest/dnntest" -T 1 -d 4 -g 1.0 -G 1.0 -H 150 -R 0.65 -i 2500 -u 48816 -x "$spec" "$R/qaplib/wil50.dat" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
export -f one; export D R
awk '{print NR, $1}' "$D/all32.txt" | xargs -P 4 -L 1 bash -c 'one "$@"' _
