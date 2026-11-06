#!/bin/sh
# Clean builds of a commit from `git archive` (no build products, no
# src/MESCHACH): with g++, with clang++, with clang for C too (DNN=0), and
# with a PATH that holds only make, g++ and what g++ itself runs.
#   git archive -o /tmp/head.tar HEAD
#   sh docs/meschach_removal/stage23/clean_build.sh /tmp/head.tar WORKDIR
# The g++ and clang++ builds make every target `make` can build from the
# repo: bnbdrv, dnntest, dnncontinuation and experiments/dnn_cuts/dnncuts.
TAR=$1; WD=$2
mkdir -p "$WD"
extract() { rm -rf "$1"; mkdir -p "$1"; tar -xf "$TAR" -C "$1"; }
build() {  # name, "extra latest/ targets", make arguments...
  d=$WD/$1; targets=$2; shift 2
  extract "$d"
  echo "== $(basename "$d"): make $*${targets:+, then make -C latest $* $targets}"
  echo "   src/ present: $( [ -e "$d/src" ] && echo yes || echo no )"
  if ( cd "$d" && make "$@" > build.log 2>&1 && { [ -z "$targets" ] || make -C latest "$@" $targets >> build.log 2>&1; } ); then
    echo "   build ok"
  else
    echo "   BUILD FAILED:"; grep -m3 -E 'error|Error' "$d/build.log" | sed 's/^/     /'; return 1
  fi
  echo "   warnings: $(grep -c 'warning:' "$d/build.log"), errors: $(grep -c 'error:' "$d/build.log")"
  echo "   compiler invocations: $(grep -oE '^(g\+\+|clang\+\+|cc|gcc|clang) ' "$d/build.log" | sort | uniq -c | tr -s ' ' | tr '\n' ',')"
  echo "   lines mentioning meschach in the build log: $(grep -ci meschach "$d/build.log")"
  for b in latest/bnbdrv latest/dnntest latest/dnncontinuation experiments/dnn_cuts/dnncuts; do
    [ -x "$d/$b" ] && echo "   built: $b"
  done
  echo "   bnbdrv links: $(ldd "$d/latest/bnbdrv" | awk '{print $1}' | tr '\n' ' ')"
}
ALL="dnntest dnncontinuation dnncuts"
build gcc "$ALL"
build clang "$ALL" CPP=clang++
build clang_cc_dnn0 "" CPP=clang++ CC=clang DNN=0
echo "== expected failure, as at afe99db: make CPP=clang++ CC=clang (DNN=1)"
d=$WD/clang_cc; extract "$d"
( cd "$d" && make CPP=clang++ CC=clang > build.log 2>&1 ) && echo "   build ok" \
  || echo "   fails: $(grep -m1 -oE 'fatal error: error in backend: Cannot select[^:]*' "$d/build.log") compiling $(grep -m1 -oE 'amx_gemm\.c' "$d/build.log")"

# a PATH with only make, g++ and the programs g++ runs (as, ld), plus sh and
# uname/rm used by the makefiles; python3 is left out, so latest/makefile
# falls back to the system LAPACK for the DNN build
bin=$WD/minpath
rm -rf "$bin"; mkdir -p "$bin"
for t in make g++ as ld sh uname rm; do ln -s "$(command -v $t)" "$bin/$t"; done
echo "== PATH with only: $(ls "$bin" | tr '\n' ' ')"
d=$WD/minpath_dnn0; extract "$d"
( cd "$d" && PATH=$bin make DNN=0 > build.log 2>&1 ) && echo "   make DNN=0: build ok, bnbdrv links $(ldd "$d/latest/bnbdrv" | awk '{print $1}' | grep -v -e vdso -e ld-linux | tr '\n' ' ')" \
  || { echo "   make DNN=0: BUILD FAILED"; tail -5 "$d/build.log"; }
d=$WD/minpath_dnn1; extract "$d"
( cd "$d" && PATH=$bin make > build.log 2>&1 ) && echo "   make (DNN=1): build ok" \
  || echo "   make (DNN=1): stops at amx_gemm.c, which needs a C compiler: $(grep -m1 -E 'not found|No such file' "$d/build.log")"
ln -s "$(command -v cc)" "$bin/cc"
( cd "$d" && PATH=$bin make > build.log 2>&1 ) && echo "   make (DNN=1) with cc added: build ok, bnbdrv links $(ldd "$d/latest/bnbdrv" | grep -oE 'lib(lapack|blas|scipy_openblas)[^ ]*' | tr '\n' ' ')" \
  || { echo "   make (DNN=1) with cc added: BUILD FAILED"; tail -5 "$d/build.log"; }
