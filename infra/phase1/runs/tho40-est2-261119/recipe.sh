# run tho40-est2-261119: written by infra/phase1/start.sh at 2026-11-19T02:19:27Z, commit 3f13899c, spec infra/phase1/specs/tho40-est.spec
# tho40, the best known 240516 (QAPSOLN/tho40.sln). INSTANCE is relative to the repository,
# PARAMS and EXTRA's paths to latest/, as in recipes/*.sh
INSTANCE=qaplib/tho40t.dat
PARAMS=../infra/phase1/recipes/tho40c.par
OPT=240516
UB=240518
EXTRA=-s\ SYM/tho40t.sym
ENV=''
COORD=--budget\ 1200\ --target-unit-cpu\ 1200\ --steal-after\ 300\ --checkpoint-every\ 60\ --ckpt-upload\ 300\ --trace
