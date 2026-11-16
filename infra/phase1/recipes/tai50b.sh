# tai50b: the whole tree, to prove or beat QAPLIB's best known 458,821,517 (Taillard's type b, n = 50;
# never proven optimal), with the cluster-region cut phase of experiments/bound_lab (README.md,
# section 1; tai50b_size/ estimates about 100 core-hours, 60-250, and 5,000-13,000 nodes). The
# settings are tai50b.par here, its regions tai50b.regions beside it. start.sh tai50b takes this
# recipe's data file and unit sizing, and the incumbent from QAPSOLN/: the best known + 1 (the
# objective values are integers, not all even: objstep.py), so the run must find 458,821,517 itself.
# Paths are relative to latest/, the solver's working directory (INSTANCE to the repository).
INSTANCE=qaplib/tai50b.dat
PARAMS=../infra/phase1/recipes/tai50b.par
OPT=458821517
UB=458821518
EXTRA=""
ENV=""
# 10-minute units, as bound_lab/laptop/tai50b.sh: a node takes 100-600 CPU s below the root; the
# root itself (about 72 minutes on one core, its phase at 6,000 iterations) runs in qapctl init
COORD="--budget 600 --target-unit-cpu 600 --steal-after 60 --checkpoint-every 120 --ckpt-upload 300 --trace"
