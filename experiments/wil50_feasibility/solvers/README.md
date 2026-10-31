# Alternatives to ADMM/PRSM for the DNN bound (research, Oct 31 2026)

Question: can a wil50 node (m = 45-48, K = (m-1)^2 + 1 = 2,117) be bounded
5-20x faster than with our PRSM, keeping the DNN bound's strength and a
certificate valid from any iterate? Three research agents, three angles;
paper hosts (arXiv and others) were blocked, so published numbers come from
search snippets and memory and are marked where unverified.

| report | angle | answer |
|---|---|---|
| `research_solvers.md` | other solvers for the same relaxation (Newton-ALM / SDPNAL+, Lagrangian-DNN with bisection or Newton bracketing as in Fujii et al., spectral and other bundle methods, boundary point, CGAL, Burer-Monteiro, restarted first-order) | No 5-20x. Every certifiable method needs a full or partial eigendecomposition per iteration; partial-spectrum methods stall because the multiplicity at the optimum is in the hundreds. Fujii et al.'s nodes cost about 10x ours. Only candidate: Halpern/restarted or Anderson-accelerated PRSM, 1.3-2.5x. |
| `research_relaxations.md` | cheaper relaxations (Ding-Wolkowicz, matrix splitting, de Klerk-Sotirov-Truetsch, grid decompositions, Lagrangian-DNN) | No. Everything cheaper is 2-10% of the objective weaker, against a 0.1% (60-unit) budget; DNN supplies QAPLIB's best lower bounds on wil50, sko49 and tho40. Dropping Y >= 0 alone costs about 2% on wil50. |
| `research_iteration.md` (+ `research_iteration_tracker.patch`) | per-iteration cost and hardware | Measured at m = 47 (530 ms/iteration on the 2.1 GHz Xeon): 56% filter GEMMs (compute-bound, ~75% of fp32 peak), 19% memory-bound transforms. CPU trims give at most 1.4x. A tight lambda_max for the Chebyshev filter (the Weyl bound is about 8x too loose) lets degree 3 beat degree 4: 5-10% fewer iterations. TF32/fp16 precision is safe for the filter, bf16 is not. The only 5-20x path is a GPU port with many nodes per device (estimated 65-130 node-iterations per device-iteration of a Xeon core, 2-8x per dollar; 3-6 weeks). |

## Bottom line

On a CPU, the remaining levers multiply to roughly 2-3x per node: an
accelerated PRSM (1.3-2.5x, untested on SDP), CPU trims (1.4x), and the tight
lambda_max (5-10%). That does not change the wil50 verdict (about 4e7 nodes;
the relaxation's strength sets the fathoming depth). A GPU port is the one
path to 5-20x per iteration and 2-8x per dollar, and only matters if a large
run is planned.

Cheapest first experiments, by report: (1) Halpern or Anderson acceleration
of PRSM behind a flag, scored on the 32 nodes of ../tuning/giveup/all32.txt;
(2) the tight lambda_max with degree 3 (`research_iteration_tracker.patch`,
DNN_TIGHT_LAM) on node 2 of all32 to convergence; (3) not needed: the
relaxations report predicts that turning off Y >= 0 alone loses ~2% at the
near misses, which already excludes everything dominated by SDP R3.
