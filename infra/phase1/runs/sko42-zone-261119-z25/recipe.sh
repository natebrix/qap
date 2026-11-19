# run sko42-zone-261119-z25: written by infra/phase1/start.sh at 2026-11-19T02:21:34Z, commit e23f534a, spec infra/phase1/specs/sko42-zone.spec
# sko42, the best known 15812 (QAPSOLN/sko42.sln). INSTANCE is relative to the repository,
# PARAMS and EXTRA's paths to latest/, as in recipes/*.sh
INSTANCE=qaplib/sko42t.dat
PARAMS=../experiments/estimates/sko42/sko42_z25.par
OPT=15812
UB=15814
EXTRA=-s\ SYM/sko42t.sym
ENV=''
COORD=--budget\ 1200\ --target-unit-cpu\ 1200\ --steal-after\ 300\ --checkpoint-every\ 60\ --ckpt-upload\ 300\ --trace
