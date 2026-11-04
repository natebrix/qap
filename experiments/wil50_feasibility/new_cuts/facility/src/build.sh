#!/bin/bash
# Build fdrv in this directory against the shared objects of latest/ (read-only).
set -e
cd "$(dirname "$0")"
L=/home/user/qap/latest
OB=$(python3 -c "import scipy_openblas32 as o; print(o.get_lib_dir())")
CF="-O3 -march=native -ffp-contract=off -DBNB_STANDALONE -DHAVE_DNN"
INC="-I$L -I/home/user/qap/src/MESCHACH"
python3 -I make_dnncut_f.py
g++ $CF $INC -c dnncut_f.C -o dnncut_f.o
g++ $CF $INC -c fdrv.C -o fdrv.o
OBJS="qpbnd.o qpbnd2.o lap.o dlap.o qap.o assign.o link1.o solution.o branch.o bnb.o util.o linalg.o fw.o estimate.o svd.o sa_qap.o stats.o console.o lgb_train.o state.o dnn.o dnncut.o amx_gemm.o"
g++ $CF fdrv.o dnncut_f.o $(for o in $OBJS; do echo $L/$o; done) -o fdrv.new -lm /home/user/qap/src/MESCHACH/meschach.a -L$OB -lscipy_openblas -Wl,-rpath,$OB
mv fdrv.new fdrv; echo built fdrv
