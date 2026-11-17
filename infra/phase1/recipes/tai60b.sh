# tai60b: the whole tree, to prove or beat QAPLIB's best known 608,215,054 (Taillard's type b, n = 60;
# never proven optimal: the best published lower bound is 592,829,061), with the cluster-region cut
# phase. experiments/bound_lab/tai60b (on claude/cool-babbage-ilnogj) certified the root at 605,315,696
# (81% of the published gap) and put the tree at 29,357 nodes (21,502-38,340) and 279 cloud-box
# core-hours (202-363; its author's range 150-600) by 284 Knuth dives from saved states. tai50b's dives
# predicted 60-250 core-hours, and its solve took about 125.
# The settings are tai60b.par here, its regions tai60b.regions beside it. start.sh tai60b takes this
# recipe's data file and unit sizing, its instance type and slots from recipes/defaults.txt (a solver
# process at n = 60 needs about 2-3 GB), and the incumbent from QAPSOLN/: the best known + 1.
# Paths are relative to latest/, the solver's working directory (INSTANCE to the repository).
INSTANCE=qaplib/tai60b.dat
PARAMS=../infra/phase1/recipes/tai60b.par
OPT=608215054
UB=608215055
EXTRA=""
ENV=""
# 10-minute units as for tai50b: a node takes several minutes near the top; the root (470 + 3,000
# iterations, hours on one core) runs in qapctl init
COORD="--budget 600 --target-unit-cpu 600 --steal-after 60 --checkpoint-every 120 --ckpt-upload 300 --trace"
