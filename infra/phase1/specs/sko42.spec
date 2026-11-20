# sko42 to the end under sko42.par as of Nov 19 (experiments/estimates/sko42/zone_nov19/README.md): cut.zone
# 0.025 (two whole-tree runs: the same trees for 0.65 of 0.010's CPU) and the bound lab's sko49 policy block
# (dnn.giveup_warm = 100, cut.handover = 27:0:5, wc=1.4 in dnn.beta_sched; experiments/bound_lab/wil50_policy),
# which sko49.par has too. For two things: what sko49 will take (the same family, the same settings), and
# the block's effect on a whole tree, paired with sko42-zone-261119-z25 (these params less the block) on the
# subtrees both close, in DNN iterations. The incumbent is start.sh's default, 15,814, as z25's. About
# 1,050-1,400 CPU h (1,400-1,700 at 0.025 without the block, less the block's 1.2-1.35x): 8-11 h on two workers.
# Nov 20: with dnn.retry = 1 (sko42.par). sko42-261120, the first run of this spec, ran without it: one failed pooled
# solve (36:7,28:0, certified -36,269 against its parent's 15,640) grew a 236,000-node, 387 CPU h subtree by 08:46 UTC.
# Without such failures the tree is about 500-700 CPU h (zone_nov19/README.md, "Failed pooled solves"): 4-5.5 h.
PROBLEM=sko42
RUN=sko42-261120b
PARAMS=experiments/estimates/sko42/sko42.par
WORKERS=2
HOURS=30
