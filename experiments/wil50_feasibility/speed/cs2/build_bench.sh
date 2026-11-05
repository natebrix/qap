#!/bin/bash
# Build rr_bench against the scipy-openblas32 wheel the solver links (latest/makefile).
D=$(cd "$(dirname "$0")" && pwd)
L=$(python3 -c "import scipy_openblas32 as o; print(o.get_lib_dir())")
g++ -O3 -march=native -ffp-contract=off "$D/rr_bench.C" -o "$D/rr_bench" \
  -L"$L" -lscipy_openblas -Wl,-rpath,"$L" -lm
