#!/bin/bash
# The Polyak reference on node 13 (21:28,8:41,28:19) with the dual test's settings
# (stronger_bounds/dual/run_five.sh, the step-7 protocol: 400 uncut iterations,
# 12 rounds of 300 with 1,000 regional triangles a round, Polyak toward best + 50
# uncapped, quiet round of 1,500 from the best round's multipliers), with the
# per-check trace (-vc) so that the certified value can be read against
# iterations. Two runs, one core each:
#   polyak_ref_13.log   the binary built before this study (OLD=path)
#   polyak_new_13.log   this tree's binary, the same flags plus -savepoolall
#                       12:pool12_13.txt: the whole pool after round 12's cuts
#                       (13,000 triangles, zero multipliers included), the
#                       frozen pool of the CP-ADMM runs. Its trace must equal
#                       the old one's (check_identity.sh).
#   bash run_ref.sh            (both; OLD defaults to the scratch build)
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../../.." && pwd)
OLD=${OLD:-/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad/dnncuts_old}
NEW=$R/experiments/dnn_cuts/dnncuts
F="-grid 5:10 -u 48816 -T 1 -k 1000 -v -I 400 -i 300 -r 12 -delta 50 -nocap -freeze 0:1500 -vc -x 21:28,8:41,28:19"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$R" || exit 1
mkdir -p "$D/logs"
one() { # tag core binary extra...
  local tag=$1 core=$2 bin=$3; shift 3
  local f=$D/logs/${tag}_13.log
  grep -q '^FINAL' "$f" 2>/dev/null && return 0
  taskset -c $core "$bin" $F "$@" qaplib/wil50.dat > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
}
case ${1:-both} in
  old) one polyak_ref ${CORE:-0} "$OLD" ;;
  new) one polyak_new ${CORE:-1} "$NEW" -savepoolall 12:$D/logs/pool12_13.txt ;;
  both) one polyak_ref 0 "$OLD" & one polyak_new 1 "$NEW" -savepoolall 12:$D/logs/pool12_13.txt; wait ;;
esac
