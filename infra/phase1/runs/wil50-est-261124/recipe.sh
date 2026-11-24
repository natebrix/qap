# run wil50-est-261124: written by infra/phase1/start.sh at 2026-11-24T16:39:43Z, commit 0a4c0f83, spec infra/phase1/specs/wil50-est.spec
# wil50, the best known 48816 (QAPSOLN/wil50.sln). INSTANCE is relative to the repository,
# PARAMS and EXTRA's paths to latest/, as in recipes/*.sh
INSTANCE=qaplib/wil50t.dat
PARAMS=../infra/phase1/recipes/wil50_d2r10_g150xl2.par
OPT=48816
UB=48818
EXTRA=-s\ SYM/wil50t.sym
ENV=''
COORD=--budget\ 1200\ --target-unit-cpu\ 1200\ --steal-after\ 300\ --checkpoint-every\ 60\ --ckpt-upload\ 300\ --trace
