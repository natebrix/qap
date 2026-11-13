# tho40, the smoke test for the tho40 recipe: one depth-3 subtree of tho40t with a
# known answer. The Nov 9 A/B (experiments/xbar_lines_sweep/tho40_ab/) solved
# 26:3,31:7,37:0 under exactly these settings (tho40.par: the nug30 recipe with
# branch.xbar_lines = 2) in 110 nodes, 107 fathomed, 14,278 FW iterations, 197 CPU s
# on the M2 Pro. Through the coordinator the tree must be the same (ship.bits = 64),
# and the CPU total prices a tho40 node on the instance's cores before the full run
# (recipes/tho40.sh). Paths are relative to latest/, the solver's working directory.
INSTANCE=qaplib/tho40t.dat
PARAMS=../experiments/estimates/tho40/tho40.par
OPT=               # a subtree: no optimum to check; the totals against the reference are the check
UB=240516          # the best known value (QAPLIB; never proven optimal)
EXTRA="-F 26:3,31:7,37:0 -s SYM/tho40t.sym"
ENV=""
COORD="--budget 60 --target-unit-cpu 60 --steal-after 20 --checkpoint-every 20 --ckpt-upload 60"
