# run sko42-est2-261121: written by infra/phase1/start.sh at 2026-11-21T17:30:35Z, commit 3044ac84, spec infra/phase1/specs/sko42-est2.spec
# sko42, the best known 15812 (QAPSOLN/sko42.sln). INSTANCE is relative to the repository,
# PARAMS and EXTRA's paths to latest/, as in recipes/*.sh
# a dives run of sko42-est-261121b (QAP_FROM, QAP_DIVES): node.sh makes it with qapctl dives; this recipe is not used
INSTANCE=qaplib/sko42t.dat
PARAMS=../experiments/estimates/sko42/sko42.par
OPT=15812
UB=''
EXTRA=-s\ SYM/sko42t.sym
ENV=''
COORD=--budget\ 1200\ --target-unit-cpu\ 1200\ --steal-after\ 300\ --checkpoint-every\ 60\ --ckpt-upload\ 300\ --trace
