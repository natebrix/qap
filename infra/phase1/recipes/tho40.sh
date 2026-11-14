# tho40: the whole tree, to prove or beat QAPLIB's best known 240,516 (Thonemann and
# Bölte, 5 x 8 grid; never proven optimal). The settings are experiments/estimates/tho40/
# (the Nov 8 nug30 record recipe with branch.xbar_lines = 2, ship.bits = 64; the Nov 8
# estimate, experiments/tho40_estimate/README.md, is 6e5 to 3e6 nodes, 0.9e6 to 4e6
# CPU s at the M2 Pro's 1.39 s per node), on tho40t.dat with the grid's symmetry (12
# depth-1 orbits). Launch with --no-reference (a bench/parallel.py rerun would double
# the work) and a --limit of days, e.g.
#   infra/phase1/launch.sh tho40-a --recipe tho40 --on-demand --type c7g.16xlarge --no-reference --limit 4320
# Paths are relative to latest/, the solver's working directory.
INSTANCE=qaplib/tho40t.dat
PARAMS=../experiments/estimates/tho40/tho40.par
OPT=240516         # the value the run is expected to end on; a lower one is a new best (compare.py then says FAIL: read it)
UB=240516          # the incumbent given to the solver
EXTRA="-s SYM/tho40t.sym"
ENV=""
# 20-minute units (about 900 nodes each at 1.4 s): with 64 slots and the default reserve
# of 2 x slots, 128 units stay queued; a unit's checkpoint goes up every 5 minutes
# --trace: every node's line (QAP_NODE_TRACE) and the DNN nodes' (QAP_DNN_TRACE), per unit in the
# store (runs/RUN/units/U/a<attempt>/nodes.txt, the frontier's in runs/RUN/frontier/), for the tree
# viewer: a few hundred MB at the estimate's 0.6e6 to 3e6 nodes, a few % of CPU on QPB nodes
COORD="--budget 1200 --target-unit-cpu 1200 --steal-after 300 --checkpoint-every 60 --ckpt-upload 300 --trace"
