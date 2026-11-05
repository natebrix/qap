#!/bin/sh
# Build a scratch bnbdrv whose my_svd can use LAPACK (MES_EIG=dsyev|dsyevd|dsyevr).
#   docs/meschach_removal/build_lapack.sh DEST
# DEST/latest and DEST/src/MESCHACH are fresh copies of this checkout's; the
# checkout itself is not modified. Needs the scipy-openblas32 wheel.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
D=$1
[ -n "$D" ] || { echo "usage: $0 DEST" >&2; exit 1; }
mkdir -p "$D/src"
rm -rf "$D/latest" "$D/src/MESCHACH"
cp -r "$ROOT/latest" "$D/latest"
cp -r "$ROOT/src/MESCHACH" "$D/src/MESCHACH"
cd "$D/latest"
rm -f ./*.o ./*.d bnbdrv
sed -i 's/^void my_svd(MAT \*A, VEC \*d)$/void my_svd_meschach(MAT *A, VEC *d)/' linalg.C
grep -q '^void my_svd_meschach' linalg.C
cat "$HERE/my_svd_lapack.C" >> linalg.C
make -C ../src/MESCHACH basic > ../meschach_build.log 2>&1
make -j2 > ../solver_build.log 2>&1
echo "built $D/latest/bnbdrv"
