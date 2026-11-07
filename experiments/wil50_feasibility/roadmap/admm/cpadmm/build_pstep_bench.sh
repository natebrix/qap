#!/bin/bash
# build pstep_bench.C against this tree's latest/ objects (make -C latest bnb first):
#   bash build_pstep_bench.sh   ->  ./pstep_bench next to this script
D=$(cd "$(dirname "$0")" && pwd)
W=$(cd "$D/../../../../../latest" && pwd)
cd "$D" || exit 1
O=""
for f in qpbnd qpbnd2 lap dlap qap assign link1 solution branch bnb util linalg fw estimate mat sa_qap stats console lgb_train state symeig dnn amx_gemm; do O="$O $W/$f.o"; done
L=/usr/local/lib/python3.13/dist-packages/scipy_openblas32/lib
g++ -O3 -march=native -ffp-contract=off -Wno-register -DBNB_STANDALONE -DHAVE_DNN -I"$W" pstep_bench.C $O -o pstep_bench -lm -L$L -lscipy_openblas -Wl,-rpath,$L
