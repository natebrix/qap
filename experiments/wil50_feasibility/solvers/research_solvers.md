# Alternative solvers for the DNN relaxation at wil50 nodes (m = 45-48): survey (Oct 31 2026)

Question: is there a solver for the *same* facially reduced DNN relaxation (K = (m-1)^2+1,
K = 2117 at m = 47) that reaches a fathoming decision 5-20x faster than our PRSM
(130 ms/iteration on M2 Pro, 300-500 ms on the 2.1 GHz Xeon; 150-700 iterations per node)?

Short answer: **no published method does, and the structure of the problem explains why.**
Every method that produces a certified bound needs, per iteration, either a full
eigendecomposition of a K x K matrix (ADMM/PRSM, APG/Lagrangian-DNN, SDPNAL+, quasi-Newton
dual ascent) or a partial one (bundle, CGAL, Burer-Monteiro). Partial-spectrum methods are
cheap per iteration but converge far too slowly for the 0.1-0.3% accuracy a fathoming
decision needs, because the DNN optimum at these nodes has rank in the hundreds (our tracker
sees hundreds of eigenvalues of each sign). So the per-iteration floor is the eigenstep (the
other agent's topic) and the only solver-side lever is the *iteration count*. The realistic
gains are 1.5-3x from acceleration/restarts of the splitting we already run. Details follow.

Web access: arxiv.org, optimization-online.org, opus4.kobv.de and ar5iv were not reachable
for fetch; figures below come from search snippets and memory, flagged where unsure.

## Per method

**1. Lagrangian-DNN + bisection / Newton-bracketing (Kim-Kojima-Toh 2016, Arima et al. 2017,
KKT OMS 2021; used by Fujii-Ito-Kim-Kojima-Shinano-Toh, ZIB 21-02 / arXiv 2101.09629).**
Relaxation: Lagrangian penalty of the equalities, cone = PSD(n^2+1) intersect (nonneg orthant
with gangster zeros); equivalent in value to the DNN. Inner solver: accelerated proximal
gradient on a projection problem, one full eigendecomposition of order n^2+1 per APG step
(no facial reduction, so N = 2210 at m = 47, slightly larger than our K). Outer: bracketing
on the scalar Lagrange multiplier; the Newton-bracketing variant needs fewer outer steps
than bisection but each outer step is hundreds to a few thousand APG iterations. Bound
validity: yes, their "robust" scheme (Arima et al. 2017) makes the lower bound valid at every
bisection iterate; the Newton-bracketing paper also gives a valid lower bound at any iterate.
Warm start: the outer bracket and the APG point can be inherited (Fujii et al. do this).
Code: BBCPOP (MATLAB), ParaQapNB not public. Speed (derived from their tables via snippets):
tai30a 3.40e7 nodes, 5.81e5 s on 1,728 cores = **about 30 core-s per node at m <= 30**;
sko42 6.02e6 nodes, 5.12e5 s on 5,184 cores = **about 440 core-s per node at m ~ 37-39**
(Xeon-class cores, 2017-2020). Scaling the eigencost by K^3 to m = 47 gives roughly
1,000-1,500 core-s per node, against our 100-130 Xeon-s (section 9 of NOTES.md). The
Newton-bracketing paper itself says it is 2-5x *slower* than BBCPOP's bisection per bound,
and BBCPOP's root bound for nug30 takes ~1,000-2,000 s. **Estimate: 0.1-0.3x of our PRSM.**
Not a candidate; our node is already about 10x cheaper than theirs. (Their wil50 node-count
estimate agrees with ours; this is why they have not solved it.)

**2. SDPNAL+ / semismooth Newton augmented Lagrangian (Zhao-Sun-Toh 2010; Yang-Sun-Toh 2015).**
Outer ALM on the dual, inner semismooth Newton-CG; each CG matvec needs products with the
eigenvector basis of the current K x K matrix, so one Newton step costs several
eigendecomposition-equivalents, and there is a long ADMM+ warm-up phase. Certification:
any outer multiplier Z plugs into our Proposition (any Z gives a bound), so yes in principle,
but intermediate Newton iterates are not better bounds than ADMM iterates at the same
eigen-budget; the method's strength is 1e-6 accuracy, which we do not need. nug30 root DNN
(matrix order 901): about 2,100 s and 2,300 iterations in BBCPOP's rerun; our C++ root
bound on nug30 is 10-20 s. Code: MATLAB only. Warm start: awkward (ALM penalty schedule).
**Estimate: well below 1x for time to a fathoming decision.** Useful only as a reference
solver for converged DNN values on a handful of nodes.

**3. Spectral bundle (Helmberg-Rendl 2000; Ding-Grimmer 2023; Liao-Ding-Zheng 2023;
Angell-McCallum 2023 "USBS", which has QAP/TSP experiments and warm starts).**
Maximizes a dual function with a lambda_max term using r extreme eigenvectors from Lanczos:
O(K^2 r) per iteration, 10-50x cheaper than a full eig at K = 2117. Certification: every
dual iterate is a valid bound (it is exactly our certificate function
s + GL(L+Z) - (m+1) lambda_max(V'ZV), maximized over Z). Problem: linear convergence needs
the bundle to capture the multiplicity of lambda_max at the optimum, which equals the rank
of R (hundreds here); with small r the method is sublinear and in practice needs thousands
of iterations for 0.1% accuracy. Rendl-Sotirov (2007) needed hours per bound at n = 30 on
2003 hardware, with interior-point inner SDPs; USBS reports QAPLIB n = 136-198 bounds whose
quality does not match SDPNAL+ accuracy. Code: ConicBundle (C++), USBS (Python/JAX).
Warm start: yes (bundle and Z). **Estimate: 0.3-1x; upside only if the negative eigenspace
at fathoming nodes is much smaller than I believe (easy to check: dnntest -v prints it).**

**4. Bundle methods for QAP SDP bounds (Rendl-Sotirov 2007).** Dualize the hard constraints
(nonnegativity, gangster) and solve the remaining small SDP by interior point in the inner
loop. Different split of the same relaxation; valid bound at every iterate; but each inner
IPM at K = 2117 is O(K^3) with a large constant, and the published runs were hours at n = 30.
**Estimate: well below 1x.** Not pursued.

**5. Boundary point / Povh-Rendl-Wiegele (2006) and OWX ADMM.** The boundary point method
is an augmented Lagrangian with one eigendecomposition per iteration, i.e. the ancestor of
the ADMM we run; same per-iteration cost, usually more iterations. **Estimate: about 1x.**

**6. Conditional gradient (CGAL, SketchyCGAL; Yurtsever et al. 2021).** One lambda_min
Lanczos per iteration (cheap), but O(1/eps) iterations and the nonnegativity constraint has
to be penalized; their QAP section uses a weaker relaxation and reports gaps of several
percent. Lower bound from an intermediate iterate: only via a dual reconstruction, loose at
their accuracy. **Estimate: not usable at the 0.1% level.**

**7. Low-rank Burer-Monteiro (LoRADS/cuLoRADS, HALLaR/cuHALLaR 2024-25, DC/rank-one
reformulations of QAP, Jiang-Zhao-Ding 2019).** Primal-side factorizations give primal
values, not lower bounds; a certified bound needs a dual reconstruction whose quality is only
guaranteed at the exact factor rank. The DNN optimum rank here is in the hundreds, so the
rank advantage mostly disappears and the Y >= 0 constraint on an N x N matrix cannot be
factorized away. **Estimate: no certified gain; not a candidate.**

**8. Accelerated / restarted first-order methods on the same splitting: Halpern-PR (HPR-LP
2024/MPC 2026, HPR-QP 2025), restarted PDHG (PDLP; PDCS 2025 has no SDP cone yet), Anderson
acceleration of the Douglas-Rachford fixed point (Ouyang et al. 2020).** Identical
per-iteration cost (one eigenstep); the multiplier Z of any iterate certifies exactly as now;
warm starts carry over unchanged. Published evidence is on LP/QP: HPR with restarts gives
O(1/k) KKT rates and typically 2-5x fewer iterations than unaccelerated PDHG/ADMM at
moderate accuracy, and Anderson-accelerated DR gives 2-3x on graphics problems. No SDP-specific
QAP numbers exist. PRSM already captured part of this (1.5x); the remaining acceleration is
not guaranteed to compound. **Estimate: 1.3-2.5x in iterations to the fathoming bound,
i.e. the same per node.** This is the only solver-side candidate worth code.

**9. Dual ascent with L-BFGS on the certificate function (BiqCrunch-style, Krislock-Malick-
Roupin 2014/2017).** Maximize s + GL(L+Z) - (m+1) lambda_max(V'ZV) directly over Z, smoothing
the polyhedral part with a quadratic penalty (their alpha) so the gradient is a projection,
which costs one full eig per function evaluation (same as an ADMM step). All iterates certify;
warm start from the parent Z is natural; BiqCrunch reports a few dozen to 150 evaluations per
node. Against our 150-700 PRSM iterations that is 1-3x, but line searches cost extra evals and
the bound of the regularized dual is weaker by the alpha term until alpha is driven to 0.
**Estimate: 1-2x, more code than (8), less certain.**

## Ranking

1. Halpern/restarted PRSM or Anderson-accelerated DR on the current splitting: 1.3-2.5x,
   cheap, certificate and warm start unchanged.
2. Spectral bundle on the certificate function with warm start: only if the measured
   dimension of the negative eigenspace of W at fathoming iterations is below ~50 (then
   O(K^2 r) per iteration could be 10x cheaper and the bundle can capture it). Check first.
3. L-BFGS dual ascent on the smoothed certificate: 1-2x, second priority after (1).

Not candidates for time to a fathoming decision: Lagrangian-DNN/Newton-bracketing (our node
is ~10x cheaper than Fujii et al.'s), SDPNAL+ (100x slower, high accuracy only), CGAL,
Burer-Monteiro, interior-point bundle.

## First experiment (candidate 1)

Prototype in `latest/dnn.C`, behind a flag: Halpern anchoring of the PRSM fixed-point pair
(Y, Z) with anchor (Y0, Z0) = warm start, step u_{k+1} = (1/(k+2)) u_0 + (1 - 1/(k+2)) T(u_k),
plus a restart rule on the relative KKT residual (restart when it drops by 0.2x, or at least
every 100 iterations; reset anchor to the current point). Alternatively, Anderson(m = 5) on
the DR fixed point with the safeguard "accept only if the certified bound does not drop".
Both are 50-100 lines; the eigenstep, certificate and GL code are untouched.
Score on the 32 nodes in `experiments/wil50_feasibility/tuning/giveup/all32.txt` with
`dnntest -T 1 -d 4 -g 1.0 -G 1.0 -H 300 -u 48816`, same machine, 4 at a time, and report
iterations and CPU to fathom (18 nodes) and the certified bound at iteration 200 on the 14
misses. Success: median iterations to fathom 290 -> <= 170 with no new false give-up and
the 200-iteration bound on the misses no worse; that is 1.7x and worth adopting. Below 1.3x,
abandon solver-side work: the remaining lever is the eigenstep or a stronger relaxation.
For candidate 2 the first step is free: `dnntest -v` on two fathoming nodes (11, 12) and two
misses, read the tracked negative-eigenspace dimension over the last 100 iterations.

Sources used (search snippets): arXiv 2101.09629 and the ParaQapNB UG-workshop slides
(tai30a/sko42 node counts, times, cores); arXiv 1905.12840 (Newton-bracketing, 2-5x slower
than BBCPOP); arXiv 1804.00761 (BBCPOP; SDPNAL+ nug30 2.1e3 s); arXiv 2006.01529 (rPRSM);
arXiv 2008.07067, 2307.07651, 2312.11801 (spectral bundle); arXiv 1912.02949 (SketchyCGAL
QAP); arXiv 2408.12179, 2507.02470 (HPR-LP/QP); Krislock-Malick-Roupin TOMS 2017; Rendl-
Sotirov Math. Prog. 2007; Arima-Kim-Kojima-Toh COAP 2017.
