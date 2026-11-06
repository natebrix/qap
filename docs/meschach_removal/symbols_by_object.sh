#!/bin/sh
# Which Meschach symbols (functions and globals defined in meschach.a) each
# solver object file references. Run from the repo root after `make` (and
# `make -C latest dnntest dnncontinuation dnncuts` for those objects).
set -e
# It needs a tree that still has Meschach (src/MESCHACH, removed from the branch
# on Nov 6): a checkout of 7b49d09, the last commit before the eigen port.
[ -f src/MESCHACH/meschach.a ] || {
  echo "$(basename "$0") needs src/MESCHACH/meschach.a: run it, after make, in a checkout of 7b49d09" >&2
  echo "  git worktree add ../qap-7b49d09 7b49d09" >&2
  exit 1; }
tmp=$(mktemp)
nm -g --defined-only src/MESCHACH/meschach.a 2>/dev/null | awk 'NF==3{print $3}' | sort -u > "$tmp"
for o in latest/*.o; do
  printf '%s: ' "$(basename "$o")"
  nm -u "$o" | awk '{print $2}' | sort -u | comm -12 - "$tmp" | tr '\n' ' '
  echo
done
rm -f "$tmp"
