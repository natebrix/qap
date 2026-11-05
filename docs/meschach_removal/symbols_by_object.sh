#!/bin/sh
# Which Meschach symbols (functions and globals defined in meschach.a) each
# solver object file references. Run from the repo root after `make` (and
# `make -C latest dnntest dnncontinuation dnncuts` for those objects).
set -e
tmp=$(mktemp)
nm -g --defined-only src/MESCHACH/meschach.a 2>/dev/null | awk 'NF==3{print $3}' | sort -u > "$tmp"
for o in latest/*.o; do
  printf '%s: ' "$(basename "$o")"
  nm -u "$o" | awk '{print $2}' | sort -u | comm -12 - "$tmp" | tr '\n' ' '
  echo
done
rm -f "$tmp"
