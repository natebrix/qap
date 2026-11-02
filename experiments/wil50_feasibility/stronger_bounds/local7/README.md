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
