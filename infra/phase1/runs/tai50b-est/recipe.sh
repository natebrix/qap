# run tai50b-est: written by infra/phase1/start.sh at 2026-11-18T19:15:38Z, commit 37d590a1
# tai50b, the best known 458821517 (QAPSOLN/tai50b.sln). INSTANCE is relative to the repository,
# PARAMS and EXTRA's paths to latest/, as in recipes/*.sh
INSTANCE=qaplib/tai50b.dat
PARAMS=../infra/phase1/recipes/tai50b.par
OPT=458821517
UB=458821518
EXTRA=''
ENV=''
COORD=--budget\ 600\ --target-unit-cpu\ 600\ --steal-after\ 60\ --checkpoint-every\ 120\ --ckpt-upload\ 300\ --trace
