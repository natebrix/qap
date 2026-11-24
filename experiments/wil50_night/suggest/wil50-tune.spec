# wil50's first tuning round: four arms side by side, each solving the same sample of wil50-est-261124's
# open nodes to the end (FROM, SAMPLE: coord/sample.py), compared node by node
# (experiments/estimates_nov19/paired.py), as tho40-tune-261119 was. The arms come from a night's screening on
# five complete subtrees (experiments/wil50_night):
#   base   wil50_d2r10.par as it is
#   g150x  cut.cp_wc_gap = 150 (a cut phase that starts within 150 of the cutoff keeps w = 0.25) and
#          cut.cp_extrap = 100:1:30:2: 0.91 of the DNN iterations on the same 413 nodes of two subtrees
#   xl2    branch.xbar_lines = 2: 631 nodes for 767 on five subtrees (two shrank, none grew; four line choices)
#   both   the two together
# The sampled nodes are at depth 3 and below, where the recipe branches by rule 9, so the arms differ in the
# bounds' cost and in rule 9's shortlist and in nothing above the sample.
# SAMPLE below is a placeholder sized from the screening (a depth-3 subtree there cost 0.2 to 1.1 CPU h): set it
# from the estimate's readout, to the depths that hold most of the rest of the tree, before committing this.
PROBLEM=wil50
RUN=wil50-tune-261125
FROM=wil50-est-261124
SAMPLE=3:60,4:60
SAMPLE_SEED=1
ARMS="base=infra/phase1/recipes/wil50_d2r10.par g150x=infra/phase1/recipes/wil50_d2r10_g150x.par xl2=infra/phase1/recipes/wil50_d2r10_xl2.par both=infra/phase1/recipes/wil50_d2r10_g150xl2.par"
VISIBILITY=private
WORKERS=1
HOURS=12
