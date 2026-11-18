# run tai60b-261117: written by infra/phase1/start.sh at 2026-11-17T05:13:43Z, commit 316361fc
# tai60b, the best known 608215054 (QAPSOLN/tai60b.sln). INSTANCE is relative to the repository,
# PARAMS and EXTRA's paths to latest/, as in recipes/*.sh
INSTANCE=qaplib/tai60b.dat
PARAMS=../infra/phase1/recipes/tai60b.par
OPT=608215054
UB=608215055
EXTRA=''
ENV=''
COORD=--budget\ 600\ --target-unit-cpu\ 600\ --steal-after\ 60\ --checkpoint-every\ 120\ --ckpt-upload\ 300\ --trace
