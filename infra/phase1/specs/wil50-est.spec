# A first estimate of wil50 under the best guess at its settings, infra/phase1/recipes/wil50_d2r10_g150xl2.par
# (Nov 24): wil50_d2r10.par (Nov 23: the bound lab's policy block, dnn.retry = 1, rule 10 at depth 2, cut.zone
# 0.012; path check PASS) with the two candidates of a night's screening on five complete subtrees, near cut
# phases at w with far ones extrapolated (0.93 of the DNN iterations) and branch.xbar_lines = 2 (631 nodes for
# 767). Nathan, Nov 24: "I'd rather go aggressive here." Neither candidate has a solved sample behind it yet,
# and under branch.xbar_lines = 2 alone the root and its 15 children keep wil50_d2r10's lines (the
# extrapolation's effect on them was not checked before this launch).
# Estimated the way tho40 and sko42 were (Knuth dives k1.5, 900 s a unit, 5 dives at least per branched open
# node, sampling by depth). wil50 has had no run under the coordinator: the Nov 7-11 sampling study sized it at
# 3,000 to 18,600 core-hours (experiments/estimates/wil50/refit.txt) under older settings.
# Four workers (sko42-261124 has finished: the Spot quota is free) and a solve phase of 1.5 h of wall (the
# workers join after about 10 minutes), about 340 CPU h: the root, the 15 depth-1 nodes and the 635 depth-2
# nodes took 71 CPU h in that study, and the 14,192-node depth-3 frontier is about 150-250 more [e], so the
# dives start from depths 3 and 4. A deep frontier on purpose: on sko42, from a 0.3 h solve, the dives read
# 0.42 of the truth from depth-2 open nodes, 0.63 from depth 3, 0.80 from depth 4, 0.60 in all (324 CPU h
# against 539; experiments/estimates_nov19, "The estimator against a known tree"). On tho40, from a 1 h solve,
# est3 read 6,844 CPU h and tho40-261122 took 6,795 (docs/next_runs.md), though under another root row. So
# read this estimate as the tree's size or somewhat under it: 1.0 to 1.7 times it.
# The decision point off: the dives run to the end whatever they cost, and the time limit bounds the run
# (an estimate cut short reads part-way, the pending units imputed). The solve stays paused in the store.
# The incumbent is start.sh's default, 48,818 (the best known 48,816 + 2).
PROBLEM=wil50
RUN=wil50-est-261124
PARAMS=infra/phase1/recipes/wil50_d2r10_g150xl2.par
ESTIMATE=1.5:900:k1.5
ESTIMATE_ENV="QAP_EST_MIN_DIVES=5 QAP_EST_SAMPLE_DEPTHS=2:1,3:0.6,4:0.2,5:0.03,6:0.03,7:0.03"
ESTIMATE_DECIDE=off
VISIBILITY=private
WORKERS=4
HOURS=10
