#!/bin/sh
# Build a timing-instrumented bnbdrv in a scratch copy of the tree.
#   docs/meschach_removal/build_timer.sh DEST
# DEST gets DEST/latest and DEST/src/MESCHACH (copies of this checkout's),
# linalg.C's my_svd gets begin/end hooks, and bnbdrv is relinked with
# mes_wrap.o and -Wl,--wrap for every Meschach symbol the solver calls.
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
cp "$HERE/mes_wrap.C" "$D/latest/"
cd "$D/latest"
rm -f ./*.o ./*.d bnbdrv
# hooks: time my_svd and optionally dump its input
sed -i 's/^void my_svd(MAT \*A, VEC \*d)$/void mes_svd_begin(MAT *); void mes_svd_end(int);\nvoid my_svd(MAT *A, VEC *d)/' linalg.C
sed -i 's/^  const unsigned int n = A->m;$/  const unsigned int n = A->m; mes_svd_begin(A);/' linalg.C
sed -i 's/^  M_FREE(H); M_FREE(Q);$/  M_FREE(H); M_FREE(Q); mes_svd_end(n);/' linalg.C
grep -q mes_svd_end linalg.C
make -C ../src/MESCHACH basic > ../meschach_build.log 2>&1
OBJS="qpbnd.o qpbnd2.o lap.o dlap.o qap.o assign.o link1.o solution.o branch.o
      bnb.o util.o linalg.o fw.o estimate.o svd.o sa_qap.o stats.o console.o
      lgb_train.o state.o dnn.o dnncut.o amx_gemm.o driver.o"
make -j2 $OBJS > ../solver_build.log 2>&1
g++ -O3 -march=native -ffp-contract=off -I../src/MESCHACH -I. -c mes_wrap.C
WRAP=""
for s in m_get m_free m_resize _m_copy m_zero m_transp v_get v_free v_resize \
         _v_copy _in_prod px_get px_free px_rows Hfactor makeHQ trieig; do
  WRAP="$WRAP -Wl,--wrap=$s"
done
OB=$(python3 -c "import scipy_openblas32 as o; print(o.get_lib_dir())")
g++ -O3 -march=native -ffp-contract=off $OBJS mes_wrap.o -o bnbdrv $WRAP -lm \
    ../src/MESCHACH/meschach.a -L"$OB" -lscipy_openblas -Wl,-rpath,"$OB"
echo "built $D/latest/bnbdrv"
