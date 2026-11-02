# The local-7 gate: is higher-order consistency violated at our iterates?

The memo's remaining candidates (the regional moment block, candidate 6, and
Astra's section-10b moment separation) assume that the regional moments at a
converged iterate lie measurably outside the higher-order consistent set. The
triple cuts (step 8) died on that assumption: violations of 0.01 where 0.25
was expected. This diagnostic asks the same question for the memo's cheap
proxy (section 3.5): for a region S and seven free facilities, the moments
(p, W) must lie in the convex hull of the admissible 0/1 configurations of
the seven membership variables (cardinality-respecting), the local marginal
polytope; the most violated inequality a.p + <B,W> <= beta with coefficients
in [-1, 1] is a 29-variable LP, on the triangle inequalities' scale (the
triangles are the 3-variable case, so the local-7 violation is at least the
triangle violation on the same facilities; the quantity of interest is the
excess).

`local7_separate.py` runs it on a dumped iterate (`dnntest -dumpY`; the
`../triples/<id>.Y.bin` files already hold the uncut 1,500-iteration iterates
of node 13 and child d2) over the 13 threshold regions, 300 random 7-subsets
plus the 7 most fractional facilities per region, and dumps the top cuts as
dense G matrices; `make_ladd7.py` assembles them at a fixed multiplier for
`dnntest -Ladd L.bin:CONST` (every cut is valid for every permutation, so the
certified value stays a bound). `run_local7.sh` does both, the trial only if
the maximum violation reaches the threshold (0.08).

Reading: maximum excess over the triangles near 0.01, as for the triples:
the regional variables are already consistent to order 7 and neither the
moment block nor 10b can gain a level here. Violations of 0.1 or more on
many subsets: run the trial, and if it gains, the higher-order route is
worth its weeks.

Plumbing check (nug12 root iterate, Nov 2, 2 s): local-7 max 0.024 on the
3-location region (a cardinality cut; no triangle is violated there), 5 cuts
dumped; through dnntest -Ladd with the cutoff at the optimum 578 they certify
568.01 at mu 1 (uncut 567.98) and 550.5 at mu 10 (overshoot), never above
578.

## Results (laptop, Nov 2; `n13.local7.{log,json}`, `d2.local7.{log,json}`, `<id>.mu<MU>.log`)

Gate passed on both nodes; the fixed-multiplier trial is positive on both.
Uncut iterates from `../triples/<id>.Y.bin` (1,500 PRSM iterations); 301
7-subsets per region (300 random + the 7 most fractional), 13 threshold
regions, 3,913 LPs in 5-6 s per node.

Per region (local-7 max violation / excess over the triangle max on the same
subsets / subsets over 0.05 / over 0.1):

| region | n13: max | excess | >0.05 | >0.1 | d2: max | excess | >0.05 | >0.1 |
|---|---|---|---|---|---|---|---|---|
| row<1 | 0.040 | 0.033 | 0 | 0 | 0.039 | 0.037 | 0 | 0 |
| row<2 | 0.116 | 0.069 | 57 | 5 | 0.117 | 0.076 | 75 | 6 |
| row<3 | 0.167 | 0.119 | 68 | 2 | 0.098 | 0.051 | 30 | 0 |
| row<4 | 0.042 | 0.031 | 0 | 0 | 0.046 | 0.023 | 0 | 0 |
| col<1 | 0.055 | 0.021 | 1 | 0 | 0.036 | 0.012 | 0 | 0 |
| col<2 | 0.034 | 0.013 | 0 | 0 | 0.041 | 0.035 | 0 | 0 |
| col<3 | 0.079 | 0.042 | 17 | 0 | 0.060 | 0.044 | 2 | 0 |
| col<4 | 0.196 | 0.136 | 91 | 11 | 0.097 | 0.057 | 40 | 0 |
| col<5 | 0.188 | 0.119 | 163 | 46 | 0.117 | 0.081 | 79 | 6 |
| col<6 | 0.146 | 0.091 | 77 | 13 | 0.129 | 0.090 | 40 | 1 |
| col<7 | 0.097 | 0.055 | 18 | 0 | 0.059 | 0.037 | 3 | 0 |
| col<8 | 0.031 | 0.015 | 0 | 0 | 0.041 | 0.026 | 0 | 0 |
| col<9 | 0.018 | 0.014 | 0 | 0 | 0.024 | 0.024 | 0 | 0 |
| all | 0.196 | 0.136 | 492 | 77 | 0.129 | 0.090 | 269 | 13 |

The violations live where the triangles' do, in the mid regions (col<4 to
col<6, row<2 and row<3, capacity 19-29 of 46-47), and they are 1.7 to 2.6
times the triangle violations on the same 7-subsets. The 40 dumped cuts
(violations 0.116-0.196 on n13, sum 5.75; 0.086-0.129 on d2, sum 3.95) all
come from those regions.

Trial, 40 cuts at fixed mu through `dnntest -Ladd`, 1,500 iterations:

| node | uncut | mu 1 | mu 5 | mu 20 | ms/iter uncut → mu 1 / 5 / 20 |
|---|---|---|---|---|---|
| 13 | 48,632.33 | 48,637.57 (+5.2) | 48,644.47 (+12.1) | 48,199.81 (−433) | 150 → 148 / 167 / 175 |
| d2 | 48,705.26 | 48,708.76 (+3.5) | 48,713.39 (+8.1) | 48,398.18 (−307) | 129 → 127 / 138 / 154 |

At mu 5 the gain is 40% of the fixed-multiplier ceiling mu x sum of
violations (28.8 and 19.8), where the triple cuts at mu 2 lost 45-49 against
a ceiling of 3.3. Per cut at fixed mu that is 0.2-0.3 units against 0.04 for
Astra's 1,000 fixed triangles at node 1 (+39.8), with the usual caveat that
40 selected cuts are not 1,000. mu 20 overshoots as every fixed multiplier
above the pool's dual does. No value near the cutoff.

Reading against the rule: the second case, not the first. The regional
variables are not consistent to order 7 at these iterates: 2-13% of random
7-subsets in the mid regions violate the local marginal polytope by more
than 0.05, the maximum excess over the triangles is 0.09-0.14, and the cuts
move the certified bound in the right direction at the right order of
magnitude. The higher-order route is worth one more step, which is the same
step the triangles took: the local-7 inequalities as a cut type in the pool
(separation is 1.4 ms per 7-subset as a 29-variable LP in scipy; a C++ port
or a heuristic over the 128 vertices is needed for in-round separation),
run through the step-7 protocol with the quiet round on nodes 13 and d2 on
top of the capped triangles. The triangles went from +39.8 fixed to +160
optimized; if the local-7 family scales the same way from its fixed +12 it
is in the +40 to +50 range at the median, which is a third of a level, not
a level. Whether that reaches d2 (18 short) is the test; node 13 (51 short)
is the harder one.

## The pool cut type (Nov 2, `latest/dnncut.C` type 6, `dnncuts -seven K[:S]`)

The local-7 inequalities are now a cut type in the library, separated inside
the rounds: for each region, S random 7-subsets of the free elements plus
S/3 subsets of the 15 most fractional ones; the moments (p_T, W_TT) are
projected onto the convex hull of the admissible configurations by away-step
Frank-Wolfe (linear convergence; plain Frank-Wolfe stalled at a residual of
0.02 and found nothing); the projection direction, scaled to max
|coefficient| 1, with the exact right-hand side over the admissible
configurations is the cut, valid for every permutation. The K most violated
enter the pool per round and take the same Polyak step and quiet round as
the triangles. On the nug12 root iterate the C++ separator finds the col<1
violation at 0.0267 against the LP's 0.0291 (the Euclidean projection normal
is not the LP's L-infinity-optimal one, a few percent weaker). Validity:
`validate_cuts.py` has a `seven` variant; the enumeration check and the
nug12 root at cutoff 578 pass. Search: `QAP_DNN_CUT_K7`, `QAP_DNN_CUT_SUB7`.

`run_pool.sh`: the step-7 protocol on nodes 13 and d2 with 300 local-7 cuts
a round on top of the 1,000 triangles (1,000 and 3,000 subsets per region),
and alone. Against the triangle-only references (58.6 and 24.2 short) the
question is whether the optimized family reaches the +40 to +50 its fixed
+12 suggests, which would fathom d2.
