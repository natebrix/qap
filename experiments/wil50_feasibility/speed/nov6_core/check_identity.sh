#!/bin/bash
# The default path is unchanged: the binary built before this study's changes
# (OLD, saved from ab937ca) and the new one, switches unset, must print the same
# trace (every "it" line and RESULT, CPU fields stripped) on nug20's root and on
# two wil50 nodes; so must DNN_CORE_REFRESH_PARTIAL=1 (bit-identical by design).
# Then the lazy certificate's own checks on the same nodes: the fathom iteration,
# the final lb and the number of Rump tests, with and without the NOSYMV bundle.
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
OLD=${OLD:?path of the old dnntest}
NEW=${NEW:-$R/latest/dnntest}
CPU=${CPU:-3}
E="OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 DNN_CHOLQR_PASSES=1 DNN_AMX=1 DNN_R32_NONCERT=1"
F="-T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40 -H 300 -i 2500 -v"
strip() { grep -E '^(it |RESULT)' "$1" | sed -E 's/ cpu [0-9.]+//; s/ cpu=[0-9.]+ ms_per_iter=[0-9.]+//'; }
cd "$R" || exit 1
mkdir -p "$D/logs/ident"
run() { # tag binary env... -- args
  local tag=$1 bin=$2; shift 2
  local envs=(); while [ "$1" != "--" ]; do envs+=("$1"); shift; done; shift
  env $E "${envs[@]}" taskset -c $CPU "$bin" $F "$@" > "$D/logs/ident/$tag.log" 2>&1
}
res() { grep ^RESULT "$1" | sed -E 's/.*lb=([^ ]+).*status=([^ ]+) iters=([0-9]+) cpu=([0-9.]+).*certs=([0-9]+).*margin=([^ ]+).*/lb \1 \2 iters \3 cpu \4 certs \5 margin \6/'; }
# name|fixes|file|extra args
for spec in "nug20|-|qaplib/nug20.dat|-i 200 -u 2570" "w19|26:27,8:38,34:28,37:21|qaplib/wil50.dat|-u 48816" "w10|39:27,30:48,25:49|qaplib/wil50.dat|-u 48816 -i 300"; do
  IFS='|' read -r name x dat extra <<< "$spec"
  xarg=(); [ "$x" != "-" ] && xarg=(-x "$x")
  run "${name}_old" "$OLD" -- "${xarg[@]}" $extra "$dat"
  run "${name}_new" "$NEW" -- "${xarg[@]}" $extra "$dat"
  run "${name}_refpart" "$NEW" DNN_CORE_REFRESH_PARTIAL=1 -- "${xarg[@]}" $extra "$dat"
  run "${name}_lazy" "$NEW" DNN_CORE_LAZYCERT=1 -- "${xarg[@]}" $extra "$dat"
  run "${name}_lazysym" "$NEW" DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 -- "${xarg[@]}" $extra "$dat"
  for v in new refpart; do
    if cmp -s <(strip "$D/logs/ident/${name}_old.log") <(strip "$D/logs/ident/${name}_$v.log"); then echo "$name: $v identical to old (trace and RESULT, CPU stripped)"; else echo "$name: $v DIFFERS from old"; diff <(strip "$D/logs/ident/${name}_old.log") <(strip "$D/logs/ident/${name}_$v.log") | head -6; fi
  done
  for v in old new refpart lazy lazysym; do echo "  $v: $(res "$D/logs/ident/${name}_$v.log")"; done
  grep -h 'lazycert exit' "$D/logs/ident/${name}_lazy.log" "$D/logs/ident/${name}_lazysym.log" | sed 's/^/  /'
done
