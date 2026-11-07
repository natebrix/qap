# Replace or wrap the PRSM for the certified DNN bound? A literature check (Nov 3, 2026)

Question: should the PRSM in `latest/dnn.C` (facially reduced DNN, K = (m-1)^2+1, about
2,117 at m = 47; 0.33 s per iteration in production; 150-1,500 iterations to fathom, give-up
near 300) be replaced or wrapped by another method, and what have others done for exactly this
relaxation? This extends `solvers/research_solvers.md` and `research_relaxations.md` and
corrects them where noted. Nothing was run. [m] = a measurement (ours, with the log or memo
named, or a published figure); [e] = my estimate. arxiv.org and math.uwaterloo.ca were not
reachable for fetch, so published figures come from search snippets of the papers' tables;
the ones I could not pin down are marked "unverified".

## Short answer

**Do not replace it.** Every published solver for this relaxation at n = 49-50 is a cold start
to 1e-5 and takes 3,400-15,300 iterations at 1.3-1.7 s each. Our warm-started, early-stopped
PRSM is 4-5x cheaper per iteration at a slightly smaller K and needs 3-30x fewer iterations
per decision (table 1). The measured spectrum rules out every partial-spectrum and low-rank
method: rank R* is about 1,500 of 2,117 and the dual has rank about 600 (`solvers/
research_iteration.md` §3, npos 1,480-1,650, nneg 460-720 [m]). Two wraps are worth a test,
both cheap to kill: **(1)** the pieces of Graham et al.'s rPRSM we never built: the
trace-constrained R-step and the doubly stochastic projection in the Y-step. `experiments/
eigensolver/README.md` lists them under "Not done". 1.0-1.5x [e], 2-3 days. **(2)** In the cut
phase only, CP-ADMM: cuts projected in the Y-step by cyclic Dykstra, as de Meijer and Sotirov
do, with the Dykstra multipliers fed to the existing certificate. Same ceiling (the dual test),
1.3-2x fewer phase iterations [e], 1-2 weeks. The GPU port stays the large lever (roadmap).

## 1. What others measured on this relaxation at n = 49-50 (root, cold start)

| source | method | wil50 lb | sko49 lb | iterations | time |
|---|---|---|---|---|---|
| Rendl-Sotirov 2007 | bundle on SDP R3, no Y >= 0 | 47,098 | 22,004 | - | hours at n = 30 [m, unverified detail] |
| Oliveira-Wolkowicz-Xu 2018 | ADMM, FR, box Y-step | 48,125 | 22,653 | 11,000-15,300 to 1e-5 (wil50, as rerun by Graham et al.) | - |
| Graham-Hu-Im-Li-Wolkowicz 2022 | rPRSM, trace on R | 48,126 | 22,654 | 4,700 (arXiv) / 5,384 (revision) | 6,133 / 9,165 s (sko49 14,235 / 13,453 s) |
| Kim-Kojima-Toh 2021 (BBCPOP and Newton-bracketing) | Lagrangian-DNN, APG, order n^2+1 | 48,122-48,123 | 22,651-22,653 | 5,453-8,650 APG (wil50); 3,383-9,170 (sko49) | 8.0e3-1.3e4 s; 4.4e3-1.2e4 s |
| ours | PRSM, tracked eigenstep, GL + Rump certificate | about 48,130 (root) | 3.1% gap | 150-1,500 per node, warm | 0.33 s/iteration |

All rows are [m]. The DNN value is the same across rows. The relaxation fixes the strength
(`research_relaxations.md`). Only the solver differs. Their seconds per iteration are 1.3-1.7
(MATLAB, a full eigendecomposition of order 2,402-2,501). Ours is 0.33 s at K = 2,117. Scaling
by K^3 makes it 3-4x on equal K and unequal hardware [e]. I found no 2023-2026 paper with new
iteration counts or bounds for this relaxation at n = 30-50. Hu-Sotirov-Wolkowicz (MP 2023)
combine FR with symmetry reduction and reach n > 500 on symmetric instances. Hou-Tang-Toh
(OR 2025, RNNAL; RiNNAL+ 2025-26) are low-rank and note that QAP DNN solutions are typically
high rank. 2403.11210 is about penalty heuristics.

## 2. Comparison at K about 2,000 (table a)

| method | cost per iteration | iterations needed | certified dual bound | memory per node | parallel / GPU | maturity |
|---|---|---|---|---|---|---|
| **PRSM, ours** | one tracked eigenstep, about 45 GF (34 GF fp32 GEMM) [m] | 150-1,500 warm [m] | Z' = beta(W - P(W)) gives V'Z'V <= 0, GL polyhedral part, Rump Cholesky every 10 iterations (143 ms, 3%) [m] | about 470 MB [m] | GEMM-bound; 3-8 ms per node-iteration on A100 [e] | production |
| rPRSM (Graham 2022) | same plus simplex projection of the eigenvalues and an O(m^2) DS projection [e] | 4,700-5,384 cold to 1e-5 [m] | dual function with box set, ceiling [m] | same | same | MATLAB research code |
| Lagrangian-DNN, APG with bisection or Newton-bracketing | one full eigendecomposition of order n^2+1 | 5,453-8,650 cold [m] | valid at every bracket step [m] | same | same | MATLAB; ParaQapNB not public |
| SDPNAL+ (SSN-ALM) | per Newton step: full eig plus 10-50 CG products at about 4K^2 min(r, K-r), about 20 GF each, i.e. 5-25 PRSM iterations [e] | ADMM+ warm-up plus tens to hundreds of Newton steps [e] | any outer Z through our certificate | 2-3x ours (Jacobian data) [e] | long dependent chains, poor batch fit [e] | MATLAB; QAP DNN is in its test set [m] |
| spectral bundle, dual (Helmberg-Rendl 2000; Ding-Grimmer 2023; USBS 2024) | Lanczos for r pairs plus a master QP with r(r+1)/2 variables | linear only if r >= mult. of lambda_max = rank R* about 1,500 [m], so a 1.1e6-variable QP; with r = 20-50, sublinear: 1e3-1e4+ [e] | every iterate is exactly our certificate function | small | Lanczos fits GPU | ConicBundle (C++), USBS (JAX) |
| spectral bundle, primal (Liao-Ding-Zheng 2023) | same | needs r >= rank of dual, about 600 [m] | primal-side | small | same | research |
| Burer-Monteiro / RNNAL / SDPLR+ / cuHALLaR | GEMMs K x K x r | r >= 1,500 [m]: no saving; Y >= 0 on 4.9M entries penalized | trace-bounded dual (SDPLR+) = our formula on a reconstructed Z | K r | good | mature for low-rank SDPs |
| CGAL / SketchyCGAL | one lambda_min | O(1/eps); the QAP runs use a 1e6-iteration cap [m] | weak; their metric is the rounding gap | tiny | good | research |
| dual-factorized ADMM (De Santis-Rendl-Wiegele 2018; Cerulli et al. 2021) | GEMMs with U, K x 600 | about ADMM's | perturbation or repair of the dual [m] | K r | good | our Chebyshev-tracked negative subspace already does this |
| RLT2 dual ascent (Hahn-Grant; Adams et al. 2007; Date-Nagi GPU) | O(n^6) multipliers, about 1.6e10 at n = 50 [e] | many | every dual-feasible LP point | far too large | GPU clusters | ">= 10x rPRSM's time for similar bounds" (Graham et al.) [m] |
| CP-ALM / CP-ADMM (de Meijer-Sotirov 2021; de Meijer et al. 2023) | PRSM plus one Dykstra sweep over the cuts, +10-30% [e] | not reported at this size | their paper needs a post-processing solve [m]; ours: GL(L + Z' + G'mu_D) - h'mu_D for the Dykstra mu_D >= 0 [e] | + pool | same as PRSM | research code, published bounds on QCCP and graph partition |
| box-cone ADMM (Charbonneau-Hu-Sotirov, arXiv 2610.41051, Oct 29 2026) | PRSM plus one sort-and-clip projection per cut row (the 2n bound-factor products of one inequality share a slack, so they form one box cone; no Dykstra loop, no slack per scalar cut) [m, their claim] | not reported on QAP; BQP against SDPNAL+, Max-Cut triangle cuts inside MADAM's B&B 6.5% median faster on 220 graphs [m] | none: "residual-corrected dual bounds ... should not be interpreted as certified" [m]; ours would be the Dykstra form of item 2 with mu_D read off the cone multipliers [e] | + pool | same as PRSM | C/MPI research code; the cheap projection exists only for cut families with the shared-slack form, which ours (pooled, not bound-factor products of a common row) may not have; an hour's read for item 2, not a build |
| preconditioned ALM / SDPF+ (Tang-Toh, arXiv 2605.17089, May 2026) | ALM with a tangent-space-weighted penalty on a feasible-method base | - | any outer Z through our certificate | - | - | claims wins only on SDPs with low-rank optimal solutions; the same exclusion as cuHALLaR, LoRADS and RNNAL (rank R* about 1,500 of 2,117) |

## 3. Why nothing replaces the PRSM

Time to a decision is cost per iteration times iterations. On cost per iteration: every method
that certifies needs the spectrum of a K x K matrix with about 1,500 positive and 600 negative
eigenvalues [m]. Methods that touch fewer eigenpairs (bundle, CGAL, Lanczos-based ALM) pay for
it in convergence. The guarantees say so: Ding-Grimmer's linear rate needs the bundle to
capture the rank of the primal solution, and Liao-Ding-Zheng's primal variant needs the rank of
the dual. Methods that touch all of them pay at least our 45 GF, and ours is already tracked,
fp32/bf16 and GEMM-bound. On iterations: no published method reports fewer than about 3,400
iterations from a cold start at n = 49-50. Our 150-1,500 comes from the warm start, the GL
certificate (1.5-6x sooner than the box bound, `dnn.C` header §2) and stopping at the cutoff.
These are properties of the use case and carry over to any splitting. A second-order method
(SDPNAL+) is the only one that could win on iterations, but at full rank on both sides each of
its Newton steps costs 5-25 of our iterations [e]. Earlier BBCPOP comparisons put SDPNAL+ at
about 2,100 s on the nug30 root DNN, against our 10-20 s (`research_solvers.md`; unverified
here).

**Corrections to the earlier surveys.**
- `research_solvers.md` item 3 left the spectral bundle open "if the negative eigenspace is
  small". That is the wrong dimension. For the dual bundle, the multiplicity of
  lambda_max(V'ZV) = 0 at the optimum equals rank R*, about 1,500 [m]. It is closed, and
  ranking item 2 needs no further check.
- `research_relaxations.md` §6 and `research_solvers.md` say "we already use PRSM" of Graham
  et al. We use the double dual step (gamma1 = gamma = 1.0; 0.9/0.9 tested [m]). We do not use
  the trace-constrained R-step or the redundant DS projection in the Y-step. Their
  "strengthened" relaxation has the same optimal value, because the trace is implied. The
  strengthening is in the bound from inexact iterates, and our (m+1) lambda_max certificate
  already has it.
- The USBS claim "QAPLIB n = 136-198, bounds below SDPNAL+ accuracy" could not be verified.
  What is verified: a 500x speedup over CGAL on a QAP SDP with about 2e9 variables. That
  relaxation is weaker than the DNN, and the speed is measured by the best rounding gap, a
  primal quantity, not a certified bound.
- The Halpern/Anderson estimate of 1.3-2.5x (item 8) was measured negative (`solvers/accel`).

## 4. Ranked recommendation (b)

**0. Keep the PRSM as the engine.** No candidate is expected to beat 0.33 s x 150-1,500
iterations. Not candidates, with the measured fact that kills each: SDPNAL+ (full rank both
sides); dual and primal spectral bundle (ranks 1,500 and 600); BM/RNNAL/SDPLR+ (rank 1,500;
the authors say QAP is high rank); CGAL family (O(1/eps)); RLT2 dual ascent (>= 10x rPRSM time,
O(n^6)); Rendl-Sotirov (R3 is 2.1% below DNN on wil50); Lagrangian-DNN (same eigenstep, cold
counts 5-9k); Kanoh-Yoshise centering ADMM (helps early iterates of cold starts, and we warm
start); smoothed L-BFGS dual ascent (BiqCrunch style: one full projection per evaluation and no
QAP evidence at this rank; 0.7-1.5x [e], last).

**1. The rPRSM pieces (wrap, 2-3 days, 1.0-1.5x [e]).** The first piece is
R = P_psd(W - tau I), with tau chosen so that trace R = m+1. tau comes from a scalar
root-find. trace(W) and the tracked eigenvalues give the sum, provided the tracked block holds
every eigenvalue below tau (the block grows by the eigenvalues in (0, tau) when tau > 0). The
second piece projects the x-row Y_0,: onto the doubly stochastic set inside the Y-step: a
Dykstra or Sinkhorn-type projection of m x m, which is negligible. The certificate stays valid
unchanged: V'Z'V = beta(tau I + (W - tau I)_-), and the Rump check bounds its lambda_max at
beta tau. The bound pays (m+1) beta tau, which goes to 0. Graham et al. measured 2.0-3.3x
fewer iterations than ADMM on wil50 [m]. We already have PRSM's 1.5x [m], so the remainder is
1.0-1.5x [e].
**Cheapest test:** the nug20 root curve (minutes), then all32 in the production configuration,
scored on median iterations to fathom and lb@200 on the misses.
**Kill:** under 10% fewer iterations, or lb@200 on the misses lower, or any fathom/miss flip.

**2. CP-ADMM for the cut phase (wrap, 1-2 weeks, 1.1-1.3x on the total [e]).** The dual test
(`stronger_bounds/dualtest`) fixes the ceiling: Polyak 12 rounds plus a quiet round is within
about 10 of the pool's dual optimum [m]. So the only gain is reaching it in fewer iterations.
The Y12s phase is about 12 x 200 + up to 1,000 quiet iterations [m]. CP-ADMM puts the frozen
pool into the Y-step: one cyclic Dykstra sweep per iteration over 5-7k cuts with 1.2-5.3k
nonzeros each, 10-30M operations, under 10% of an iteration [e]. It certifies with
GL(L + Z' + G'mu_D) - h'mu_D, valid for any mu_D >= 0 by the LP duality of
`research_bound_quality.md` §3.3. If it reaches the plateau in 1,000-2,000 iterations, the
phase is 1.3-2x faster, and at its 46% share that is 1.1-1.3x on the total [e]. The risk is
the closed fast60 proxy: frequent multiplier moves were worse. Dykstra is the convergent
version of that idea, but the heavily overlapping regional cuts may make the sweep slow to mix.
**Cheapest test:** node 13's round-12 pool, frozen. PRSM plus the Dykstra Y-step from the
give-up point.
**Kill:** not within 10 of 48,763 by 2,000 iterations, or more than 30% per-iteration overhead.

**3. Already on the roadmap and consistent with the literature:** certificate polishing / GL
subgradient ascent on mu (a bundle on the polyhedral part only, which is where Rendl-Sotirov's
idea still applies), and the GPU port (all splittings are GEMM-bound; none is better suited to
a GPU than PRSM).

## 5. What the QAP-DNN-ADMM papers did that we have not (c)

| paper | item | status here |
|---|---|---|
| OWX 2018 | FR, box/gangster Y-step, dual-projection bound | done; our bound is stronger (GL, trace, Rump) |
| OWX 2018 | low-rank/rank-one upper bounds from Y | not done; not needed (incumbent 48,816) |
| OWX 2018 | cutting planes added to the Y-subproblem ("trivial with FR") | not done; we dualize cuts in rounds. Item 2 |
| OWX 2018 | high-accuracy runs (1e-5, up to 40k iterations) | deliberately not |
| Graham 2022 | PRSM double dual step | done (1.0/1.0; 0.9/0.9 tested) |
| Graham 2022 | trace constraint in the R-step | built and measured Nov 5 (`../../speed/rprsm`): -1% iterations on wil50; its multiplier is a flat direction of our GL certificate |
| Graham 2022 | DS projection in the Y-step | not in their Algorithm 3.1, only their Remark 2.9; built Nov 5 (x-row and arrow): -10% iterations cold, -6.5% warm, -4 to -6% CPU (`../../speed/rprsm`) |
| Graham 2022 | "restricted" dual update | checked against the arXiv PDF: both dual steps projected onto Z_0 (diagonal and arrow fixed at -L); built Nov 5, neutral on nug20 (`../../speed/rprsm`) |
| Graham 2022 | ceiling of the bound for integer data | done (even objective, cutoff 48,814.05) |
| Hu-Sotirov-Wolkowicz 2023 | symmetry reduction after FR | n/a below the root: wil50's order-4 grid group fixes almost no location, so a depth-3 node's stabilizer is trivial or order 2 [e] |
| Kim-Kojima-Toh / Fujii 2021 | warm-started brackets in B&B, valid bound at any iterate | equivalent (warm Z, any-Z certificate) |
| de Meijer-Sotirov 2021/2023 | cuts in the projection by Dykstra | **not done**. Item 2 |
| Charbonneau-Hu-Sotirov 2026 | bound-factor RLT cuts as one box cone per inequality, projected by sort and clip inside the ADMM | not done; found Nov 6 by a citation pass over Graham 2022 (OpenAlex and Semantic Scholar, 14 citing papers, all others off-topic or already here); the question for item 2 is whether any of our cut families has the shared-slack form |
| Jansson-Chaykin-Keil 2007; VSDP; SDPLR+ 2024 | rigorous bound from a dual point plus a primal eigenvalue/trace bound | done, with the exact trace m+1 (stronger than a lambda_max(X) bound) and Rump's Cholesky; the certificate is 3% of time with margin about 1e-7 [m], so nothing cheaper or stronger is left on this side |

## 6. References (d)

- Rendl, Sotirov. Bounds for the QAP using the bundle method. Math. Prog. 109 (2007) 505-524.
- Helmberg, Rendl. A spectral bundle method for SDP. SIAM J. Optim. 10 (2000) 673-696; ConicBundle.
- Oliveira, Wolkowicz, Xu. ADMM for the SDP relaxation of the QAP. Math. Prog. Comp. 10 (2018) 631-658; arXiv 1512.05448.
- Graham, Hu, Im, Li, Wolkowicz. A restricted dual Peaceman-Rachford splitting method for a strengthened DNN relaxation for QAP. INFORMS J. Comput. 34 (2022) 2125-2143; arXiv 2006.01529.
- Hu, Sotirov, Wolkowicz. Facial reduction for symmetry reduced semidefinite and doubly nonnegative programs. Math. Prog. 200 (2023) 475-529.
- Kim, Kojima, Toh. A Newton-bracketing method for a simple conic optimization problem. Optim. Methods Softw. 36 (2021); arXiv 1905.12840 (wil50/sko49/tho40 rows).
- Fujii, Ito, Kim, Kojima, Shinano, Toh. Solving challenging large scale QAPs. ZIB 21-02, arXiv 2101.09629.
- Yang, Sun, Toh. SDPNAL+. Math. Prog. Comp. 7 (2015); arXiv 1406.0942.
- Ding, Grimmer. Revisiting spectral bundle methods. SIAM J. Optim. 33 (2023) 1305-1332.
- Liao, Ding, Zheng. Spectral bundle methods for primal and dual SDPs. arXiv 2307.07651 (2023); Liao, Zheng, bundle-based ALM, arXiv 2502.08835 (2025).
- Angell, McCallum. Fast, scalable, warm-start SDP with spectral bundling and sketching (USBS). ICML 2024; arXiv 2312.11801.
- Yurtsever, Tropp, Fercoq, Udell, Cevher. Scalable semidefinite programming. SIMODS 2021; arXiv 1912.02949.
- Hou, Tang, Toh. Low-rank ALM for DNN relaxations of MBQPs (RNNAL). Oper. Res. 2025; arXiv 2502.13849. RiNNAL+, arXiv 2507.13776.
- Huang, Gleich. Suboptimality bounds for trace-bounded SDPs (SDPLR+). arXiv 2406.10407 (2024).
- Charbonneau, Hu, Sotirov. Box-cone representations of bound-factor RLT constraints in semidefinite relaxations. arXiv 2610.41051 (Oct 2026).
- Tang, Toh. A preconditioned augmented Lagrangian method for solving semidefinite programming problems. arXiv 2605.17089 (May 2026).
- De Santis, Rendl, Wiegele. Using a factored dual in ALMs for SDP. Oper. Res. Lett. 2018; Cerulli, De Santis, Gaar, Wiegele. 4OR 19 (2021) 415-448.
- de Meijer, Sotirov. SDP-based bounds for the quadratic cycle cover problem via cutting plane ALM. INFORMS J. Comput. 2021; arXiv 2009.03616. de Meijer, Sotirov, Wiegele, Zhao. Partitioning through projections. Comput. Oper. Res. 151 (2023).
- Kanoh, Yoshise. Centering ADMM for the SDP relaxation of the QAP. arXiv 2001.05739 (2020).
- Hahn, Zhu, Guignard, Hightower, Saltzman. A level-3 RLT bound for the QAP. INFORMS J. Comput. 2012; Adams, Guignard, Hahn, Hightower. EJOR 180 (2007).
- Jansson, Chaykin, Keil. Rigorous error bounds for the optimal value in SDP. SIAM J. Numer. Anal. 46 (2007) 180-200; VSDP.
