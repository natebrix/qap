# run tho40-tune-261119-cross5: written by infra/phase1/start.sh at 2026-11-19T19:22:31Z, commit f076d8f2, spec infra/phase1/specs/tho40-tune.spec
# tho40, the best known 240516 (QAPSOLN/tho40.sln). INSTANCE is relative to the repository,
# PARAMS and EXTRA's paths to latest/, as in recipes/*.sh
INSTANCE=qaplib/tho40t.dat
PARAMS=../infra/phase1/recipes/tho40c_cross5.par
OPT=240516
UB=240518
EXTRA=-s\ SYM/tho40t.sym
ENV=DNN_ADMM_WARM_BETA=0.1:30:0.36:1
COORD=--budget\ 1200\ --target-unit-cpu\ 1200\ --steal-after\ 300\ --checkpoint-every\ 60\ --ckpt-upload\ 300\ --trace
