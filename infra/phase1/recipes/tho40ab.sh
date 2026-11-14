# The parent-state A/B (Nov 14), before tho40 runs whole again: one depth-2 subtree of tho40t,
# 26:10 then 8:0 (the depth-2 sweep's first-solve gap 3.00%, about 3,200 nodes by the estimate's
# fit, ~2.5 CPU-hours on c7a), solved twice with tho40's settings and 5-minute units, so that
# many nodes move between units: tho40ab-ship0 (no parent DNN state in state files: a node
# handed to another unit starts cold) and tho40ab-ship32 (fp32, one copy per parent in the
# run's store, coord/warm.py). Compare DNN iterations and CPU per node, wall, and the disk and
# store the second uses. Not launched by itself: the two recipes below set PARAMS.
#   infra/phase1/launch.sh tho40ab-0 --recipe tho40ab-ship0 --on-demand --type c7a.16xlarge --no-reference --limit 120 --disk 100
#   infra/phase1/launch.sh tho40ab-32 --recipe tho40ab-ship32 --on-demand --type c7a.16xlarge --no-reference --limit 120 --disk 100
# Paths are relative to latest/, the solver's working directory.
INSTANCE=qaplib/tho40t.dat
PARAMS=
OPT=               # a subtree: no optimum to check
UB=240516          # the best known value: both arms prune alike
EXTRA="-F 26:10,8:0 -s SYM/tho40t.sym"
ENV=""
COORD="--budget 300 --target-unit-cpu 300 --steal-after 60 --checkpoint-every 60 --ckpt-upload 120 --trace"
