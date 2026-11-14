# tho40: the whole tree, to prove or beat QAPLIB's best known 240,516 (Thonemann and
# Bölte, 5 x 8 grid; never proven optimal). The settings are tho40.par here:
# experiments/estimates/tho40/tho40.par (the Nov 8 nug30 record recipe with
# branch.xbar_lines = 2; the Nov 8 estimate, experiments/tho40_estimate/README.md, is 6e5
# to 3e6 nodes, 0.9e6 to 4e6 CPU s at the M2 Pro's 1.39 s per node) with ship.bits = 32
# through the run's parent-state store and the DNN bound at every depth (dnn.depths =
# 1:38), on tho40t.dat with the grid's symmetry (12 depth-1 orbits). tho40-a (Nov 14,
# ship.bits = 64 inline) filled its 32 GB disk at 06:31 and was stopped. Launch with
# --no-reference (a bench/parallel.py rerun would double the work), a --limit of days and
# --disk 100 (infra/phase1/README.md, Disk and failures), e.g.
#   infra/phase1/launch.sh tho40-b --recipe tho40 --on-demand --type c7a.16xlarge --no-reference --limit 4320 --disk 100
# c7a (AMD, one worker per core, 64 on a 16xlarge): nug30 at 16 workers took 0.69 CPU s a node
# on a c7a.4xlarge (runs/nug30-c7a), 1.35 on a c7g.4xlarge (runs/nug30-a), 0.34 on the M2 Pro.
# Paths are relative to latest/, the solver's working directory.
INSTANCE=qaplib/tho40t.dat
PARAMS=../infra/phase1/recipes/tho40.par
OPT=240516         # the value the run is expected to end on; a lower one is a new best (compare.py then says FAIL: read it)
# The incumbent given to the solver: the best known + 2, the next value on tho40's objective
# grid (its values are even), as the 2000 nug30 solve was given the optimum + 2. The run
# must then find a solution of 240,516 itself (status.json: found_by_search) and prove that
# none is better. The extra tree, nodes whose bound lies within 2 of the best known, is
# about 0.1% by the estimate's gap fit (nug12 with -u 580: 254 nodes against 241).
UB=240518
EXTRA="-s SYM/tho40t.sym"
ENV=""
# 20-minute units (about 900 nodes each at 1.4 s): with 64 slots and the default reserve
# of 2 x slots, 128 units stay queued; a unit's checkpoint goes up every 5 minutes
# --trace: every node's line (QAP_NODE_TRACE) and the DNN nodes' (QAP_DNN_TRACE), per unit in the
# store (runs/RUN/units/U/a<attempt>/nodes.txt, the frontier's in runs/RUN/frontier/), for the tree
# viewer: a few hundred MB at the estimate's 0.6e6 to 3e6 nodes, a few % of CPU on QPB nodes
COORD="--budget 1200 --target-unit-cpu 1200 --steal-after 300 --checkpoint-every 60 --ckpt-upload 300 --trace"
