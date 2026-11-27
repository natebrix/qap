# The first run on worker groups (Nov 19): nug30 at its optimum, one Spot worker. About 1,300 CPU s
# of work, so 20-40 minutes in all with the coordinator's build and gates; a dollar or two.
# What to look for (collect.sh nug30-smoke-261119): the coordinator's steps say "worker group
# qap-RUN: 1 instances" after init and "0 instances" at the end; the worker's steps say "left the
# worker group"; the group's activities show one launch and one termination.
PROBLEM=nug30
RUN=nug30-smoke-261119
PARAMS=infra/phase1/recipes/nug30.par
INCUMBENT=6124
FLAGS=none
UNITS="--budget 60 --target-unit-cpu 60 --steal-after 20 --checkpoint-every 20 --ckpt-upload 60"
VISIBILITY=private
WORKERS=1
HOURS=2
DISK=32
