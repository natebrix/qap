#!/bin/bash
# Shared by the instance settings files here (INSTANCE/settings.sh), which set
# FULLRUN_INSTANCE, FULLRUN_SYM (may be empty), FULLRUN_UB, FULLRUN_PAR and
# FULLRUN_GOLDEN and then source this. The names are the full run's
# (../wil50_feasibility/fullrun/settings.sh) so that the sampling scripts, which
# read a settings file through fullrun_check (sampling/common.py, SAMPLING_SETTINGS),
# take any of these files unchanged.
#
#   fullrun_dump BINARY   the binary's resolved settings under FULLRUN_PAR (the
#                         settings/setting lines of --print-settings, less par.file)
#   fullrun_check [BIN]   true when the files exist and the dump equals FULLRUN_GOLDEN
# Regenerate a golden only in a commit that means to change that instance's settings:
#   (. experiments/estimates/INSTANCE/settings.sh; fullrun_dump latest/bnbdrv > "$FULLRUN_GOLDEN")

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

fullrun_dump() {
  "$1" --print-settings -p "$FULLRUN_PAR" 2>/dev/null | grep -E '^settings? ' | grep -v '^setting par\.file '
}
fullrun_check() {
  local b=${1:-$FULLRUN_ROOT/latest/bnbdrv} ok=0
  [ -f "$FULLRUN_INSTANCE" ] || { echo "settings: $FULLRUN_INSTANCE is missing"; ok=1; }
  [ -z "$FULLRUN_SYM" ] || [ -f "$FULLRUN_SYM" ] || { echo "settings: $FULLRUN_SYM is missing"; ok=1; }
  [ -f "$FULLRUN_PAR" ] || { echo "settings: $FULLRUN_PAR is missing"; return 1; }
  [ -x "$b" ] || { echo "settings: $b is missing: make solver"; return 1; }
  [ -f "$FULLRUN_GOLDEN" ] || { echo "settings: $FULLRUN_GOLDEN is missing"; return 1; }
  if ! diff <(fullrun_dump "$b") "$FULLRUN_GOLDEN" >/dev/null; then
    echo "settings: $b's settings under $FULLRUN_PAR differ from $FULLRUN_GOLDEN:"
    diff <(fullrun_dump "$b") "$FULLRUN_GOLDEN" | sed 's/^/  /'
    ok=1
  fi
  return $ok
}
