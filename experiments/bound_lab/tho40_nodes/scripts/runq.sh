#!/bin/bash
# runq.sh JOBFILE WORKERS: each line of JOBFILE is "TAG<TAB>command..."; runs the command from the
# repository root with stdout+stderr to logs/TAG.log (skipped when logs/TAG.done exists).
D=$(cd "$(dirname "$0")/.." && pwd)
ROOT=$(cd "$D/../../.." && pwd)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 D ROOT
run1() {
  tag="${1%%$'\t'*}"; cmd="${1#*$'\t'}"
  [ -e "$D/logs/$tag.done" ] && exit 0
  cd "$ROOT" && echo "# $cmd" > "$D/logs/$tag.log" && bash -c "$cmd" >> "$D/logs/$tag.log" 2>&1
  echo $? > "$D/logs/$tag.done"
}
export -f run1
grep -v '^#' "$1" | grep . | tr '\n' '\0' | xargs -0 -P "${2:-1}" -I{} bash -c 'run1 "$1"' _ {}
