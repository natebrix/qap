#!/bin/bash
# runseq.sh JOBFILE [PID]: wait for PID to exit (if given), then run the jobs of JOBFILE one after
# another ("TAG<TAB>command", from the repository root, output to logs/TAG.log; skipped when
# logs/TAG.done exists). The file is re-read before every job, so lines can be appended while it runs.
D=$(cd "$(dirname "$0")/.." && pwd)
ROOT=$(cd "$D/../../.." && pwd)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
[ -n "$2" ] && tail --pid="$2" -f /dev/null
while :; do
  next=""
  while IFS= read -r line; do
    case "$line" in ''|'#'*) continue;; esac
    tag="${line%%$'\t'*}"
    [ -e "$D/logs/$tag.done" ] || [ -e "$D/logs/$tag.log" ] && continue
    next="$line"; break
  done < "$1"
  [ -z "$next" ] && break
  tag="${next%%$'\t'*}"; cmd="${next#*$'\t'}"
  ( cd "$ROOT" && echo "# $cmd" > "$D/logs/$tag.log" && bash -c "$cmd" >> "$D/logs/$tag.log" 2>&1; echo $? > "$D/logs/$tag.done" )
done
