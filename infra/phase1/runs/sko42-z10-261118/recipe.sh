# run sko42-z10-261118: written by infra/phase1/start.sh at 2026-11-18T21:28:02Z, commit bbf15b89
# sko42, the best known 15812 (QAPSOLN/sko42.sln). INSTANCE is relative to the repository,
# PARAMS and EXTRA's paths to latest/, as in recipes/*.sh
INSTANCE=qaplib/sko42t.dat
PARAMS=../experiments/estimates/sko42/sko42.par
OPT=15812
UB=15814
EXTRA=-s\ SYM/sko42t.sym
ENV=''
COORD=--budget\ 600\ --target-unit-cpu\ 600\ --steal-after\ 60\ --checkpoint-every\ 120
