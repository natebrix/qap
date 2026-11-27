# run nug30-smoke-261119: written by infra/phase1/start.sh at 2026-11-19T01:42:19Z, commit 1017e6cd, spec infra/phase1/specs/nug30-smoke.spec
# nug30, the best known 6124 (QAPSOLN/nug30.sln). INSTANCE is relative to the repository,
# PARAMS and EXTRA's paths to latest/, as in recipes/*.sh
INSTANCE=qaplib/nug30.dat
PARAMS=../infra/phase1/recipes/nug30.par
OPT=6124
UB=6124
EXTRA=''
ENV=''
COORD=--budget\ 60\ --target-unit-cpu\ 60\ --steal-after\ 20\ --checkpoint-every\ 20\ --ckpt-upload\ 60
