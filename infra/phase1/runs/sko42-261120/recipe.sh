# run sko42-261120: written by infra/phase1/start.sh at 2026-11-20T01:06:30Z, commit 0855ec91, spec infra/phase1/specs/sko42.spec
# sko42, the best known 15812 (QAPSOLN/sko42.sln). INSTANCE is relative to the repository,
# PARAMS and EXTRA's paths to latest/, as in recipes/*.sh
INSTANCE=qaplib/sko42t.dat
PARAMS=../experiments/estimates/sko42/sko42.par
OPT=15812
UB=15814
EXTRA=-s\ SYM/sko42t.sym
ENV=''
COORD=--budget\ 1200\ --target-unit-cpu\ 1200\ --steal-after\ 300\ --checkpoint-every\ 60\ --ckpt-upload\ 300\ --trace
