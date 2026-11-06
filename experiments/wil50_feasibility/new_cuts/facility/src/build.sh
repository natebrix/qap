#!/bin/bash
# Build fdrv in this directory: fdrv.C plus dnncut_f.C (latest/dnncut.C with
# its symbols renamed by make_dnncut_f.py), linked with the objects of a solver
# tree that has been built with make and make dnntest (read-only).
#
#   src/build.sh [LATEST]     LATEST: that tree's latest/ (default: this checkout's)
#
# make_dnncut_f.py patches dnncut.C as it was at 7eabced (Nov 4; the copy is
# dnncut_f.C.orig), and dnncut.C has changed since, so build against that commit:
#   git worktree add ../qap-7eabced 7eabced
#   (cd ../qap-7eabced && make && make dnntest)
#   experiments/wil50_feasibility/new_cuts/facility/src/build.sh ../qap-7eabced/latest
# The objects, include paths and libraries are read from LATEST's makefile, so
# the script builds against a tree with Meschach (before Nov 6) or without it.
set -e
L=$(cd "${1:-$(dirname "$0")/../../../../../latest}" && pwd)
cd "$(dirname "$0")"
H=$(pwd)
if ! cmp -s "$L/dnncut.C" dnncut_f.C.orig; then
  echo "$L/dnncut.C is not the dnncut.C fdrv was written against (dnncut_f.C.orig, from 7eabced):"
  echo "build against a checkout of 7eabced, as the top of $0 says"
  exit 1
fi
mk() {   # the value of a variable of LATEST's makefile
  make -s --no-print-directory -C "$L" -f makefile -f - _print V="$1" << 'EOF'
_print: ; @echo $($(V))
EOF
}
CF="$(mk CF) -DBNB_STANDALONE -DHAVE_DNN"
INC="$(mk INCLUDES) -I$L"
OBJS="$(mk OBJS) $(mk DNN_OBJ)"
for o in $OBJS; do [ -f "$L/$o" ] || { echo "$L/$o is missing: run make and make dnntest in that tree first"; exit 1; }; done
python3 -I make_dnncut_f.py
# compile and link from LATEST, where its makefile's relative paths (Meschach's,
# before Nov 6) point
cd "$L"
g++ $CF $INC -c "$H/dnncut_f.C" -o "$H/dnncut_f.o"
g++ $CF $INC -c "$H/fdrv.C" -o "$H/fdrv.o"
g++ $CF "$H/fdrv.o" "$H/dnncut_f.o" $OBJS -o "$H/fdrv.new" $(mk LDLIBS) $(mk DNN_LAPACK)
mv "$H/fdrv.new" "$H/fdrv"; echo "built $H/fdrv against $L"
