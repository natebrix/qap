#!/bin/sh
# Differential test of latest/mat.C against Meschach (mat_difftest.C), run
# from the repo root:
#   sh docs/meschach_removal/stage23/mat_difftest.sh OLD [STEPS] [WORK]
# OLD is a source tree that still has src/MESCHACH with meschach.a built (for
# example a `git archive afe99db` extracted and built with make). WORK
# (default: a temporary directory) gets the binaries and outputs.
set -e
OLD=$1; STEPS=${2:-400000}; WORK=${3:-$(mktemp -d)}
T=docs/meschach_removal/stage23/mat_difftest.C
MES=$OLD/src/MESCHACH
g++ -O2 -DUSE_MESCHACH -I"$MES" "$T" "$MES/meschach.a" -lm -o "$WORK/t_mes" 2>/dev/null
g++ -O3 -march=native -ffp-contract=off -Ilatest "$T" latest/mat.C -o "$WORK/t_mat"
echo "mat.C: $(g++ --version | head -1), -O3 -march=native -ffp-contract=off; Meschach: g++ -O2 test driver, meschach.a from $MES"
for seed in 1 2 3 4 5 6 7 8; do
  "$WORK/t_mes" $seed "$STEPS" > "$WORK/a.txt"
  "$WORK/t_mat" $seed "$STEPS" > "$WORK/b.txt"
  ops=$(grep -c ' op ' "$WORK/a.txt"); skips=$(grep -c '^skip' "$WORK/a.txt" || true)
  if cmp -s "$WORK/a.txt" "$WORK/b.txt"; then r=identical; else r=DIFFERENT; fi
  echo "seed $seed: $ops operations ($skips skipped), outputs $r, $(wc -c < "$WORK/a.txt") bytes, md5 $(md5sum < "$WORK/a.txt" | cut -c1-12) / $(md5sum < "$WORK/b.txt" | cut -c1-12)"
done
echo "mat.C built by clang++ -O3 -march=native -ffp-contract=off ($(clang++ --version | head -1)):"
clang++ -O3 -march=native -ffp-contract=off -Ilatest "$T" latest/mat.C -o "$WORK/t_matc"
for seed in 1 2 3; do
  "$WORK/t_mes" $seed "$STEPS" > "$WORK/a.txt"; "$WORK/t_matc" $seed "$STEPS" > "$WORK/b.txt"
  if cmp -s "$WORK/a.txt" "$WORK/b.txt"; then echo "  seed $seed: identical"; else echo "  seed $seed: DIFFERENT"; fi
done
# valgrind does not decode AVX-512, so its mat.C build is plain -O2
g++ -O2 -Ilatest "$T" latest/mat.C -o "$WORK/t_mat_vg"
for b in t_mes t_mat_vg; do
  echo "valgrind --leak-check=full $b 1 20000:"
  valgrind --leak-check=full "$WORK/$b" 1 20000 2>&1 > "$WORK/$b.out" | grep -E 'Invalid|definitely lost:|indirectly lost:|ERROR SUMMARY' | sed 's/^==[0-9]*==/ /'
done
cmp -s "$WORK/t_mes.out" "$WORK/t_mat_vg.out" && echo "  (outputs identical; the lost blocks are VECs made by v_get(0), which v_resize and v_copy replace rather than resize, in both)"
