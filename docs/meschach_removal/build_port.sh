#!/bin/sh
# Build a scratch bnbdrv whose my_svd uses sym_eig_port (eig_port.C), the
# plain-array port of Meschach's Hfactor/makeHQ/trieig, instead of Meschach.
#   docs/meschach_removal/build_port.sh DEST
# The checkout itself is not modified.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
D=$1
[ -n "$D" ] || { echo "usage: $0 DEST" >&2; exit 1; }
mkdir -p "$D/src"
rm -rf "$D/latest" "$D/src/MESCHACH"
cp -r "$ROOT/latest" "$D/latest"
cp -r "$ROOT/src/MESCHACH" "$D/src/MESCHACH"
cp "$HERE/eig_port.C" "$HERE/eig_port.h" "$D/latest/"
cd "$D/latest"
rm -f ./*.o ./*.d bnbdrv
sed -i 's/^void my_svd(MAT \*A, VEC \*d)$/void my_svd_meschach(MAT *A, VEC *d)/' linalg.C
grep -q '^void my_svd_meschach' linalg.C
cat >> linalg.C <<'EOF'

#include "eig_port.C"
void my_svd(MAT *A, VEC *d)
{
  d = v_resize(d, A->m);
  sym_eig_port(A->m, A->me, d->ve);
}
EOF
make -C ../src/MESCHACH basic > ../meschach_build.log 2>&1
make -j2 > ../solver_build.log 2>&1
echo "built $D/latest/bnbdrv"
