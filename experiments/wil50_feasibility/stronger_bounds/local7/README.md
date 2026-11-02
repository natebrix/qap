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

## Pool results (laptop, Nov 2; `<proto>_<id>.log`, two jobs at a time, one core each)

Negative on both nodes: the local-7 cuts in the pool, with or without the
triangles, end 15-25 below the triangle-only references after the quiet
round. No fathom on d2; node 13 loses ground.

| node | protocol | best Polyak round | after quiet | vs triangles | drops | round cpu (s) | total cpu (s) |
|---|---|---|---|---|---|---|---|
| 13 | triangles only (ref, `../dual/tri_13.log`) | 48,708.97 (r12) | 48,755.42 | — | 2 | 62-66 | 1,067 |
| 13 | triseven (1,000 tri + 300 local-7, 1k subsets) | 48,693.63 (r10) | 48,738.16 | −17.3 | 3 | 45-50 | 866 |
| 13 | triseven3k (3k subsets) | 48,688.64 (r12) | 48,732.67 | −22.8 | 2 | 46-51 | 872 |
| 13 | sevenonly (300 local-7, no triangles) | 48,686.70 (r11) | 48,730.49 | −24.9 | 2 | 47-51 | 879 |
| d2 | triangles only (ref, `../dual/tri_d2.log`) | 48,730.62 (r8) | 48,789.87 | — | 2 | 53-62 | 1,021 |
| d2 | triseven | 48,715.42 (r12) | 48,771.95 | −17.9 | 2 | 40-44 | 752 |
| d2 | triseven3k | 48,718.28 (r12) | 48,775.26 | −14.6 | 2 | 39-44 | 755 |
| d2 | sevenonly | 48,716.33 (r12) | 48,772.63 | −17.2 | 3 | 39-45 | 766 |

Every run is 400 + 12 x 300 Polyak iterations and a 1,500-iteration quiet
round from the best round's multipliers; the quiet round adds 44-57 in every
case, as in the references. Nothing came within 38 of the cutoff; no
validity question arises.

The cuts do bite. On node 13 the separator finds 10,859 violated 7-subsets of
about 17,300 at round 0 (max 0.166; 31,948 of 52,000 with 3k subsets) and
4,000-4,800 at round 12 (max 0.06-0.08); on d2 9,736 (max 0.122) falling to
4,300-4,700 (max 0.05-0.08). The local-7 maximum falls from 0.17 to 0.06
while the triangle maximum falls from 0.09 to 0.04-0.06, so the local-7
violations keep their 1.5-2x ratio to the triangles' throughout; the 3,000
subsets per region find larger maxima at the same K = 300 but add nothing.

The 300 local-7 cuts a round alone (sevenonly) are worth 84-87% of 1,000
triangles a round: +134/+116 against +159/+133 above the uncut value. Per
cut they are 2.8x the triangles, which matches the fixed-multiplier ratio
(0.2-0.3 a cut against 0.04). But the two families together are no better
than either alone, and worse than the triangles: the best Polyak round is
15 lower in every joint run, and the quiet round, which adds the same 45-55
on top, does not recover it. The joint pool is a superset of the triangle
pool, so its dual optimum is at least the triangles', and the dual test
(step 9) put the Polyak protocol within about 10 of the triangle pool's
optimum; the loss is therefore in the dynamics, not in the family. The
likely mechanism is the one the fixed-mu trial showed at mu 20: a local-7
cut moves 35 moments at coefficient up to 1 where a triangle moves 6, so
the same multiplier step perturbs the PRSM iterate several times as much,
and 300 inner iterations recover less of it (the fixed-mu trial already
needed 1,500 iterations to realize 40% of the ceiling at mu 5). The
violations the local-7 cuts express are largely the triangles' violations
seen through a larger window (the gate found them in the same regions, on
the same elements), so the two families spend the Polyak step on the same
slack.

Reading: the pool cut type does not reach the +40 to +50 the fixed +12
suggested, and the higher-order route as a Lagrangian family does not add
to the triangles under the step-7 protocol. Two cheap checks remain before
closing it: (a) the local-7 cuts at a small fixed multiplier on top of the
triangle pool's best multipliers (additivity at fixed mu, 1,500 iterations,
one run per node), and (b) the joint pool under the slow30 protocol of step
9, which approximates the dual optimum regardless of step dynamics; if (b)
is also at or below the triangle pool's slow30 value (48,762 / 48,795) the
family is closed. A smaller K (50-100 local-7 cuts a round) or adding them
only after round 6 (`-addto`) are variants of the same question and not
worth separate runs unless (a) is positive.
