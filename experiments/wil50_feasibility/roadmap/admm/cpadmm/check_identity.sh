#!/bin/bash
# Opt-in only: with the CP-ADMM switch off (DNNParams::ycuts NULL, no -cpadmm), the
# binaries built from this tree must print the same traces as the ones built
# before the study (scratch builds of the parent commit; OLD_DNNTEST, OLD_DNNCUTS).
# As speed/nov6_core/check_identity.sh: every "it" line and RESULT, CPU stripped.
#   1. dnntest, plain settings and production settings (DNN_AMX=1
#      DNN_CHOLQR_PASSES=1 DNN_R32_NONCERT=1, -d 3 -L 2 -g 1 -G 1 -r 40, plus the
#      production bundle DNN_BETA_SCHED/LAZYCERT/NOSYMV/NOSYMW/REFRESH_PARTIAL), on
#      nug20's root, wil50 node 19 (26:27,8:38,34:28,37:21) and node 10 (39:27,30:48,25:49).
#   2. the cut phase: dnncuts with the dual test's flags on node 13, old against
#      new, from run_ref.sh's two logs (the new one also wrote the round-12 pool):
#      every "it", "round" and "quiet" line, CPU stripped.
#   3. the CP code path with an empty pool against the plain solve from the same
#      give-up point (dnncuts -cpadmm with pool "-" against -cpplain): the
#      restructured Y/Z step reproduces the plain one bit for bit.
# One core (CPU, default 3). Logs in logs/ident/.
D=$(cd "$(dirname "$0")" && pwd); R=$(cd "$D/../../../../.." && pwd)
S=/tmp/claude-0/-home-user-qap/f096d970-f4ab-52d1-8c70-e0605bd8cf8a/scratchpad
OLD_DNNTEST=${OLD_DNNTEST:-$S/dnntest_old}
NEW_DNNTEST=${NEW_DNNTEST:-$R/latest/dnntest}
NEW_DNNCUTS=${NEW_DNNCUTS:-$R/experiments/dnn_cuts/dnncuts}
CPU=${CPU:-3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$R" || exit 1
L=$D/logs/ident; mkdir -p "$L"
strip() { grep -E '^(it |RESULT|round|quiet)' "$1" | sed -E 's/ cpu [0-9.]+//; s/ cpu=[0-9.]+ ms_per_iter=[0-9.]+//'; }
same() { # label fileA fileB
  if cmp -s <(strip "$2") <(strip "$3"); then echo "$1: identical ($(strip "$2" | wc -l) lines)"; else echo "$1: DIFFERS"; diff <(strip "$2") <(strip "$3") | head -6; fi
}
PROD_ENV="DNN_AMX=1 DNN_CHOLQR_PASSES=1 DNN_R32_NONCERT=1"
BUNDLE="DNN_BETA_SCHED=two,f=1.6667,t=60,w=0.25,wc=0.7 DNN_CORE_LAZYCERT=1 DNN_CORE_NOSYMV=1 DNN_NOSYMW=1 DNN_CORE_REFRESH_PARTIAL=1"
PF="-T 1 -d 3 -L 2 -g 1.0 -G 1.0 -r 40 -H 300 -v"
for spec in "nug20|-|qaplib/nug20.dat|-i 200 -u 2570" "w19|26:27,8:38,34:28,37:21|qaplib/wil50.dat|-u 48816 -i 300" "w10|39:27,30:48,25:49|qaplib/wil50.dat|-u 48816 -i 300"; do
  IFS='|' read -r name x dat extra <<< "$spec"
  xarg=(); [ "$x" != "-" ] && xarg=(-x "$x")
  for arm in plain prod bundle; do
    case $arm in
      plain) envs=""; flags="-T 1 -v" ;;
      prod) envs="$PROD_ENV"; flags="$PF" ;;
      bundle) envs="$PROD_ENV $BUNDLE"; flags="$PF" ;;
    esac
    for v in old new; do
      bin=$OLD_DNNTEST; [ $v = new ] && bin=$NEW_DNNTEST
      f=$L/${name}_${arm}_$v.log
      [ -s "$f" ] || env $envs taskset -c $CPU "$bin" $flags "${xarg[@]}" $extra "$dat" > "$f" 2>&1
    done
    same "dnntest $name $arm" "$L/${name}_${arm}_old.log" "$L/${name}_${arm}_new.log"
  done
done
if [ -s "$D/logs/polyak_ref_13.log" ] && [ -s "$D/logs/polyak_new_13.log" ]; then
  same "dnncuts node 13 cut phase (12 rounds + quiet)" "$D/logs/polyak_ref_13.log" "$D/logs/polyak_new_13.log"
else echo "dnncuts node 13 cut phase: run_ref.sh logs missing"; fi
for spec in "nug12s|-grid 3:4 -c 578 -I 200|$D/logs/validate/nug12_swapped.dat|400" "w13|-grid 5:10 -u 48816 -I 40 -x 21:28,8:41,28:19|qaplib/wil50.dat|60"; do
  IFS='|' read -r name flags dat its <<< "$spec"
  [ -s "$dat" ] || { echo "$name: $dat missing (run validate_cp.py first)"; continue; }
  echo "-" > "$L/empty.pool"
  [ -s "$L/${name}_cpplain.log" ] || taskset -c $CPU "$NEW_DNNCUTS" $flags -T 1 -v -cpplain -cpiters $its "$dat" > "$L/${name}_cpplain.log" 2>&1
  [ -s "$L/${name}_cpempty.log" ] || taskset -c $CPU "$NEW_DNNCUTS" $flags -T 1 -v -cpadmm "$L/empty.pool" -cpiters $its "$dat" > "$L/${name}_cpempty.log" 2>&1
  same "CP path, empty pool, vs plain: $name" "$L/${name}_cpplain.log" "$L/${name}_cpempty.log"
done
