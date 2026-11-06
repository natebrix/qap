#!/bin/sh
# AddressSanitizer + UndefinedBehaviorSanitizer (leak checking on) on the quick
# suite and a brute-force subset, for two commits built the same way, so that
# only reports that are new can be singled out. Meschach (old commit) is
# instrumented too, as mat.C is in the new one. Run from the repo root:
#   git archive -o /tmp/old.tar afe99db; git archive -o /tmp/new.tar HEAD
#   sh docs/meschach_removal/stage23/sanitize.sh /tmp/old.tar /tmp/new.tar WORKDIR
# ASan's reports (leaks at exit included) go to files under WORKDIR/logs_*
# (with -fsanitize=address,undefined it is UBSAN_OPTIONS's log_path that
# counts). UBSan's runtime errors go to the solver's stderr whatever log_path
# says, so they are looked for in bench.py's per-instance logs (stdout +
# stderr), and for brute force UBSan is made fatal (halt_on_error), which
# brute_force.py counts as a crash, as it does a leak or ASan error (nonzero
# exit status).
OLD=$1; NEW=$2; WD=$3
SAN="-fsanitize=address,undefined -fno-omit-frame-pointer -g"
CF="-O1 -march=native -ffp-contract=off $SAN"
REPO=$(pwd)
for side in old new; do
  d=$WD/$side
  if [ ! -x "$d/latest/bnbdrv" ]; then
    rm -rf "$d"; mkdir -p "$d"
    if [ $side = old ]; then tar -xf "$OLD" -C "$d"; else tar -xf "$NEW" -C "$d"; fi
    (
      cd "$d"
      if [ -d src/MESCHACH ]; then
        make -C src/MESCHACH basic CFLAGS="-O -std=gnu17 -ffp-contract=off $SAN" > build.log 2>&1
      fi
      # -Wno-register: harmless for g++, already in the new makefile
      make CF="$CF -Wno-register" >> build.log 2>&1
    ) || { echo "$side: build failed"; tail "$d/build.log"; exit 1; }
  fi
  echo "$side: built with CF=\"$CF\"$( [ -d "$d/src/MESCHACH" ] && echo ", Meschach with the same sanitizers")"
done

summarize() {  # log directory
  n=$(ls "$1" | wc -l)
  echo "   ASan report files: $n"
  [ "$n" -eq 0 ] && return
  cat "$1"/* | grep -E '^SUMMARY' | sort | uniq -c | sed 's/^/   /'
  # leaks: first frame outside the sanitizer runtime and the allocator
  cat "$1"/* | awk '/^(Direct|Indirect) leak/{want=1; next} want && /#[0-9]+ / && !/libasan|asan_|malloc|calloc|realloc|operator new/{sub(/^ *#[0-9]+ 0x[0-9a-f]+ in /,""); print "   leak from: " $0; want=0}' | sed -E 's#/[^ ]*/(latest|src)/#\1/#' | sort | uniq -c
}

for side in old new; do
  d=$WD/$side
  logs=$WD/logs_$side; rm -rf "$logs"; mkdir -p "$logs/quick" "$logs/bf"
  echo "== $side: quick suite (bench.py check --tree)"
  ASAN_OPTIONS="detect_leaks=1:log_path=$logs/quick/asan" UBSAN_OPTIONS="print_stacktrace=1:log_path=$logs/quick/asan" \
    OPENBLAS_NUM_THREADS=1 python3 "$REPO/bench/bench.py" check --tree "$d" 2>&1 | tail -1 | sed 's/^/   /'
  summarize "$logs/quick"
  last=$(ls -td "$REPO"/bench/results/check-* | head -1)
  echo "   solver logs with a UBSan/ASan message ($(ls "$last/logs" | wc -l) logs in bench/results/$(basename "$last")): $(grep -l -E 'runtime error|Sanitizer' "$last"/logs/* | wc -l)"
  echo "== $side: brute force --sizes 4-6 --count 2 --jobs 2, UBSan fatal"
  ASAN_OPTIONS="detect_leaks=1:log_path=$logs/bf/asan" UBSAN_OPTIONS="print_stacktrace=1:halt_on_error=1:log_path=$logs/bf/asan" \
    OPENBLAS_NUM_THREADS=1 python3 "$REPO/tests/brute_force.py" --solver "$d/latest/bnbdrv" --sizes 4-6 --count 2 --jobs 2 2>&1 | tail -1 | sed 's/^/   /'
  summarize "$logs/bf"
done
