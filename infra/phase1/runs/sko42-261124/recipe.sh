# run sko42-261124: written by infra/phase1/start.sh at 2026-11-24T02:32:43Z, commit f801ef08
# sko42, the best known 15812 (QAPSOLN/sko42.sln). INSTANCE is relative to the repository,
# PARAMS and EXTRA's paths to latest/, as in recipes/*.sh
INSTANCE=qaplib/sko42t.dat
PARAMS=../infra/phase1/recipes/sko42_d2r10_sb.par
OPT=15812
UB=15814
EXTRA=-s\ SYM/sko42t.sym
ENV=''
COORD=--budget\ 1200\ --target-unit-cpu\ 1200\ --steal-after\ 300\ --checkpoint-every\ 60\ --ckpt-upload\ 300\ --trace
