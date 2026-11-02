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

## Results (laptop, Accelerate, one core each, Nov 2; logs `n13.*.log`, `d2.*.log`)

Negative at every multiplier on both nodes. Cutoff 48,814.05.

| node | m | uncut 1,500 it. | triples violated | violation range (sum) | mu 2 | mu 10 | mu 50 | tri + best mu, after quiet | tri alone (dual/) |
|---|---|---|---|---|---|---|---|---|---|
| 13 (`21:28,8:41,28:19`) | 47 | 48,632.33 | 60/60 (46 used) | 0.0002..0.162 (1.45) | 48,587.22 (−45.1) | 47,772.08 (−860) | 42,314.51 | 48,700.25 (113.8 short) | 48,755.42 (58.6 short) |
| d2 (`39:27,30:48,25:49,26:5`) | 46 | 48,705.26 | 60/60 (55 used) | 0.0004..0.113 (1.64) | 48,656.77 (−48.5) | 47,670.45 (−1,035) | 41,125.66 | 48,736.20 (77.9 short) | 48,789.87 (24.2 short) |

Separation: 4.0 s (d2) and 6.1 s (n13) per triple on average, max 9.6 and
14.0 s. Per iteration the dense term costs +5-10% at mu 2 and +20-25% at mu
50 in dnntest (n13: 150 -> 159 -> 184 ms/iter); the dnncuts rounds with the
term ran at 47-63 s against 60-75 s for the triangle-only runs, which were
done with four concurrent jobs rather than two, so that comparison says
nothing. No value approached the cutoff; validity is not in question.

Reading. The trial's design bounds what it could have shown: with fixed
multipliers the gain over the uncut relaxation is at most mu times the sum of
the violations at the uncut optimum, which is 3.3 (n13) / 3.3 (d2) at mu 2,
16 at mu 10 and 82 at mu 50. The observed losses of 45-50 at mu 2 are an
order of magnitude larger than the available gain, so the penalty's dense
term (580-670k nonzeros) sets back the 1,500-iteration certificate by more
than any triple cut could add; at mu 10 and 50 the overshoot is the same
effect Astra saw at mu 3. On top of the triangles the term costs 54-55 units
after the quiet round on both nodes.

The structural problem is the violation scale. Random triples at these
converged iterates violate by 0.0002-0.16 (median about 0.01), not Astra's
0.25-0.28 (ten triples at the 1,500-iteration node-1 iterate; whether those
were selected or random is not recorded). Even at her scale, a +130 level from
this family would need a sum of violations around 130/mu, i.e. hundreds of
selected triples at mu near 1 (separation hours per node at 4-14 s each), or
a multiplier that the inner solve cannot absorb. A bundle on the multipliers
would recover the overshoot but not the scale. Candidate 5 is retired for
wil50 unless a selection rule finds triples with violations an order of
magnitude larger than random ones; the 60 random violations per node are in
the `.sep.log` files for anyone who wants to look for such a rule.
