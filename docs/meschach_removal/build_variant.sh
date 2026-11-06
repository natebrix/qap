#!/bin/sh
# Build bnbdrv in a scratch copy with Meschach compiled by another compiler
# or with other flags; the solver itself is built as usual.
#   docs/meschach_removal/build_variant.sh DEST CC "CFLAGS"
# e.g. build_variant.sh /tmp/fma clang "-O -march=native"   (FMA contraction on)
#      build_variant.sh /tmp/gcc3 gcc "-O3 -march=native -ffp-contract=off"
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
# It needs a tree that still has Meschach (src/MESCHACH, removed from the branch
# on Nov 6): a checkout of 7b49d09, the last commit before the eigen port.
[ -f "$ROOT/src/MESCHACH/makefile" ] || {
  echo "$(basename "$0") needs a tree with Meschach: run it from a checkout of 7b49d09" >&2
  echo "  git worktree add ../qap-7b49d09 7b49d09; ../qap-7b49d09/docs/meschach_removal/$(basename "$0") ..." >&2
  exit 1; }
D=$1; MCC=$2; MCFLAGS=$3
[ -n "$D" ] && [ -n "$MCC" ] || { echo "usage: $0 DEST CC CFLAGS" >&2; exit 1; }
mkdir -p "$D/src"
rm -rf "$D/latest" "$D/src/MESCHACH"
cp -r "$ROOT/latest" "$D/latest"
cp -r "$ROOT/src/MESCHACH" "$D/src/MESCHACH"
rm -f "$D"/latest/*.o "$D"/latest/*.d "$D"/latest/bnbdrv "$D"/src/MESCHACH/*.o "$D"/src/MESCHACH/meschach.a
make -C "$D/src/MESCHACH" basic CC="$MCC" CFLAGS="$MCFLAGS" > "$D/meschach_build.log" 2>&1
make -C "$D/latest" -j2 > "$D/solver_build.log" 2>&1
echo "built $D/latest/bnbdrv (Meschach: $MCC $MCFLAGS)"
