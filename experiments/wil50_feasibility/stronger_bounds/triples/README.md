# Triple-consistency (projected RLT2) cuts: the fixed-multiplier trial

Memo candidate 5 (research_bound_quality.md): for three free facilities the
pairwise blocks of Y must be the marginals of one joint placement; Astra's
separator (`triplet_separate.py`, from the branch research/regional-dnn-bounds,
copied with an added `--n` option; `region_cuts.py` is its helper) finds a
violated projected RLT2 inequality <G, Y> <= 0 by an LP, rounds the potentials
to multiples of 1/1024 and checks them on every injective placement, so the
cut is valid independently of the LP solver. Astra found 10 of 10 triples
violated by 0.25-0.28 at wil50 node 1 but never used a cut in a bound.

The cheap trial (no new cut type in the library): dump the DNN iterate
(`dnntest -dumpY`, `dnncuts -dumpY`), separate K random triples, add
`mu * sum G_k` as an external objective term (`make_ladd.py` ->
`dnntest -Ladd L.bin` or `dnncuts -Ladd L.bin` on top of the regional
triangles), read the certified bound. Every value printed stays a valid lower
bound: the term goes through the same `L_add` hook as the regional cuts, its
constants (0 here) through the shift, and its magnitude into the rounding
margin. Fixed multipliers measure a lower bound on the family's value (the
triangles gave +39.8 fixed and +121 optimized on node 1).

`run_trial.sh <id> <spec> [K] [n]`: uncut 1,500-iteration reference with the
iterate dumped, K triples separated, mu in {2, 10, 50} with dnntest, then the
step-7 dnncuts protocol with the best mu on top of the triangles. One core,
about 70 minutes at m = 47.

Plumbing check (nug12 root, cutoff 578 = the optimum, Nov 2): four triples
separated at the 200-iteration iterate, violations 0.008-0.213 (0.01-0.5 s
each); uncut certified 567.90, with the four cuts at mu 1: 568.05 (+0.14, the
same order as the regional cuts' +0.13 at this root), at mu 5: 556.6 (a fixed
multiplier that large overshoots, as Astra's mu 3 did); dnncuts with the term
loads it and runs. No value above the optimum.

To run on the laptop after `git pull` and `make -C latest dnntest dnncuts`
(needs scipy):
  bash experiments/wil50_feasibility/stronger_bounds/triples/run_trial.sh n13 21:28,8:41,28:19
  bash experiments/wil50_feasibility/stronger_bounds/triples/run_trial.sh d2 39:27,30:48,25:49,26:5
