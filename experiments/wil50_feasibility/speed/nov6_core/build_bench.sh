#!/bin/bash
# Build the microbenchmarks of this study against the solver's objects and BLAS
# (scipy-openblas32, as latest/makefile links it). Run from anywhere; binaries land here.
#   test_amx: ../../amx/test_amx_gemm.c against the current latest/amx_gemm.o (runtime knobs
#             DNN_CORE_AMX_IB/KB/NT/ORDER are read by the kernel itself)
#   bench_refresh: the full refresh, ssyevd against ssytrd + sstedc + a partial sormtr
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/../../../.." && pwd)
OB=$(python3 -c "import scipy_openblas32 as o; print(o.get_lib_dir())")
set -e
cd "$R/latest" && nice -n 19 taskset -c 0 make amx_gemm.o > /dev/null
cd "$D"
nice -n 19 taskset -c 0 gcc -O3 -march=native -o test_amx "$R/experiments/wil50_feasibility/amx/test_amx_gemm.c" "$R/latest/amx_gemm.o" -L"$OB" -lscipy_openblas -Wl,-rpath,"$OB" -lm
nice -n 19 taskset -c 0 gcc -O3 -march=native -o bench_refresh bench_refresh.c -L"$OB" -lscipy_openblas -Wl,-rpath,"$OB" -lm
echo built: test_amx bench_refresh
