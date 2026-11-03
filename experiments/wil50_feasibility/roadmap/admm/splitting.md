# The DNN PRSM against the general ADMM literature: what is left to take

Research memo, Nov 3 2026: which ideas from the general ADMM literature (Boyd et al. 2011
and its descendants) would cut iterations, or the cost of an iteration, of our certified
DNN bound. Nothing was run. **[m]** = measured
(file or paper given), **[e]** = my estimate. stanford.edu, arXiv and PMLR were unreachable
from this session: papers were checked through search snippets, and the SCS and OSQP rules
were read from source (SCS 3.2.7 `src/scs.c`; OSQP 0.6.7 `auxil.c`, from PyPI). Claims whose
numbers could not be retrieved are marked qualitative.

## The iteration, in our notation

Per iteration (`latest/dnn.C`, main loop):

1. W = V'(Y + Z/β)V; R = P_psd(W) by the tracked subspace iteration (inexact; W − R ⪯ 0 exactly).
2. Certificate multiplier Zc = Z + β(Y − VRV'), so V'Zc V = β(W − R) ⪯ 0 (checked every 10).
3. Z ← Z + γ1 β(Y − VRV').
4. Y ← Π_P(VRV' − (L/scale + Z)/β), P = {0 ≤ Y ≤ 1, Y00 = 1, gangster zeros}, entrywise.
5. Z ← Z + γ β(Y − VRV').

Production: γ = γ1 = 1, β = 0.3·m/3 (4.7 at m = 47) fixed; children start from the
parent's mapped Z with **Y reset to zero**.

## (a) Diagnosis

On the axes that general solvers treat as free wins, our PRSM is already at or past the
state of the art. Equilibration (Ruiz in OSQP, SCS, COSMO) is vacuous: R ↦ VRV' is an
isometry (V'V = I), and the R-step stays an eigendecomposition only under a scalar
Frobenius penalty. Over-relaxation is at its endpoint: (γ1, γ) = (1, 1) is Peaceman-Rachford
on the dual, beyond the α = 1.6 of OSQP and COSMO and the 1.5 of SCS, and on the boundary of
the He-Ma-Yuan (2016) domain. Approximate warm-started projections (Rontsis et al. 2022) are
what the tracker does. Chordal decomposition (CDCS, COSMO) has nothing to decompose (rank R
about 1,500 of 2,117), and the uncut problem has two blocks. What every production conic
ADMM has and ours lacks is **state management**: (i) a residual-adapted penalty (OSQP,
SCS 3, COSMO, SDPAD, the boundary-point method), where ours is fixed and hand-tuned,
although β off by 3.3x costs 2-3x in iterations (nug20/nug30, `dnn_bound/cpp_report.md`)
[m], the optimum differs between nug30 (0.21) and wil50 (0.3) [m], and kra30 does not
converge at the default [m]; (ii) a **full** primal-dual warm start (OSQP and SCS warm-start
x, y and s; we reset Y to zero, which costs 17-47 certified units at 300 iterations on the
same node [m]); (iii) R-step inexactness that is constant in time, where theory asks for
loose-early and tight-late (one-pass bf16 tracks fp32 to iteration 100, then drifts [m]).
The residuals agree: on wil50 node 2 (`solvers/accel/node2_d3L2.log`) pres/(β·dY), with dY
bounding the dual residual's norm from above, is 0.32-0.45 at iterations 10-50, 0.27 at 100,
0.15 at 150 and 0.118-0.125 from 200 to 1,500. The transient and the tail want different
balances, and a fixed β serves only the tail. Anderson, Halpern and momentum are closed for
diagnosable reasons (end of (b)).

## (b) Ranked ideas

| # | idea | expected effect | days | risk | evidence | cheapest first test, and the result that kills it |
|---|---|---|---|---|---|---|
| 1 | **Full-state child warm start**: map the parent's final Y to the child as Z is mapped already, instead of Y = 0 | 10-25% fewer iterations on every warm-started child [e]. Same node: restoring Y with Z gives +17 to +47 certified at 300 iterations over Z alone [m], about 50-70 iterations of node 2's trajectory at 230-300 [e] | 1-1.5 | low: Y is only the starting iterate | `same_node_continuation/RESULTS.md` [m]; OSQP warm start: 2.9-3.6x fewer iterations on portfolio back-tests [m] (Stellato et al. 2020, Table 7); warm-started B&B children (Stellato et al. 2018, qualitative) | Parents 9 and 6 of `tuning/giveup/all32.txt` run to their give-up with `save`, then their 16 depth-4 children (nodes 17-32) with Z only against Y+Z (variants A and B in (c)), `-H 0 -i 600`. **Kill:** median iterations to fathom down less than 10% **and** lb@100 up less than 10 in both variants |
| 2 | **Normalized residual-balancing β**, changed only at full refreshes and frozen after iteration 200 | wil50: 0-15% fewer iterations [e] (β is already tuned there). Instances where β is mis-scaled: up to the 2-3x that hand tuning bought on nug [m] | 2 | medium: perturbs the tracker, and the right target ratio is unknown (Wohlberg 2017) | OSQP: ρ ← ρ·sqrt(r̂/ŝ), adopted at 5x (code); SCS 3: sqrt of the geometric mean of r̂/ŝ, ≥ 100 iterations apart, adopted outside [1/√10, √10] (code); He-Yang-Wang 2000; Boyd et al. 2011 §3.4.1; Lorenz and Tran-Dinh 2019 (summable increments) | 10 lines: two-phase β on all32, 0.5 then 0.3 (×m/3) from the first refresh after iteration 60 (the predicted direction: the early ratio is 3x the tail's), and 0.2 then 0.3 as control. **Kill on wil50:** neither moves lb@200 on the misses by +20 nor the fathom median by 10 iterations. Keep for robustness only if the controller lifts the kra30a root bound above QPB (68,600; fixed β: 17,368) |
| 3 | **Scheduled R-step inexactness**: one-pass bf16 filter products while the residual is large, the 3-pass split from then on | −13% CPU per iteration in the early phase, 5-7% per node [e] | 1-2 | low: the certificate does not depend on the eigensolver | bf16 vs fp32: +0.4 at it 50, +0.9 at 100, −3.4 at 150, −7 at 200 [m] (`research_iteration.md` §2); one pass 7.4-7.8 ms vs 17-20 ms for the 3-pass kernel [m] (`amx/README.md`); Rontsis et al. 2022: summable-error approximate projections, up to 20x faster, iterations not noticeably higher [m]; Eckstein and Yao 2018 | node 2 to 300 iterations and then all32, switching at pres < 0.25·pres_10 or iteration 120. **Kill:** a lb@200 deficit above 2 on any miss, any later fathom, or less than 4% CPU per node |
| 4 | **Homogenizing scale t** (Ỹ = DYD, D = diag(t, 1, …, 1)): the one diagonal preconditioner that keeps R a projection | 0-20% [e], unknown | 3 | medium: touches forward/expand (V_D = orth(DV) changes only V's column through coordinate 0) and the certificate's trace term (trace Ỹ = t² + m) | equilibration in OSQP, SCS and COSMO: qualitative, LP/QP only, never on SDP | Zero code first: split pres² and dY² into row/column 0 and the rest at the checks. **Kill:** row/column 0 carries under 5% of both after iteration 50 |
| 5 | **Spectral (Barzilai-Borwein) penalty**, AADMM | as #2, higher variance [e] | 2 on top of #2 | medium-high: curvature estimates of indicators are degenerate | Xu, Figueiredo, Goldstein 2017: fewer iterations than vanilla, fast and residual-balanced ADMM, incl. DIMACS SDPs (qualitative; table not retrieved) | Only if #2 survives: the same all32 run. **Kill:** no gain over #2 |

**Closed, with the reason** (not re-proposed):

- **Over-relaxation and the (γ1, γ) pair.** The He-Ma-Yuan (2016) domain is γ1 ∈ (−1, 1),
  γ ∈ (0, (1+√5)/2), γ1 + γ > 0, |γ1| < 1 + γ − γ². Over it, γ1 + γ < 1 + 2γ − γ² ≤ 2, with
  the supremum at (1, 1): production already takes the largest total dual step the domain
  allows. (0.5, 1.618) lies outside (1 + φ − φ² = 0) and diverged in our tests [m]. The
  strictly contractive (0.9, 0.9) (He, Liu, Wang, Yuan 2014) lost [m].
- **Ruiz equilibration, per-constraint ρ** (OSQP: 10³ρ on equalities): our equalities sit
  in the exact Y-projection, and a non-scalar metric on the coupling term makes the R-step
  a weighted projection with no eigen-solution. Only #4 survives.
- **Chordal or multi-block splitting** (CDCS, COSMO; sGS-ADMM): dense W, two blocks; only
  dualized cuts need it (`stronger_bounds/research_bound_quality.md` §2.1).
- **Anderson, Halpern, momentum, PDLP-style restarted averaging.** PRSM here contracts
  geometrically (gap ratio 0.34-0.55 per 50 iterations on fathoms, `tuning/README.md`), so
  an anchor with O(1/k) worst-case rate can only slow it, as measured (1.9x); momentum
  duplicates PR's own extrapolation [m]; Anderson needs an exact fixed-point map, and with
  exact ssyevd it gained only 0.1 [m]. No variant is worth a run.

## (c) The top three in detail

### 1. Full-state child warm start

Today `dnn_warm_child` substitutes x_ij = 1 into Z and `dnn_bound` starts from Y = e0e0'.
With Z near the parent's multiplier, V'ZV is nearly negative semidefinite, so the first
R-step sees W ≈ V'ZV/β ⪯ 0, returns R ≈ 0, and dozens of iterations rebuild the positive
part (rank about 1,500) the parent already had. The parent's final Y is already saved
(`save->Y`), so no extra memory is needed.

Exact change. Parent P (size m), child C fixing (i, j), r = 1 + i·m + j, child index
a ↔ parent pa (the map in `dnn_warm_child`), x_r = Y_P[r][r]:

- **Variant A (restriction):** Y_C[0][0] = 1; Y_C[0][a] = Y_C[a][0] = Y_P[0][pa];
  Y_C[a][b] = Y_P[pa][pb].
- **Variant B (conditioned):** if x_r ≥ 0.05, q_a = clip(Y_P[r][pa]/x_r, 0, 1), otherwise
  q_a = Y_P[0][pa]. Then Y_C[0][a] = Y_C[a][a] = q_a, and for a ≠ b,
  Y_C[a][b] = clip(Y_P[pa][pb]·sqrt(q_a q_b / (Y_P[pa][pa] Y_P[pb][pb])), 0, 1). This is a
  congruence scaling, so diag(Y) = x holds as it does on the DNN feasible set.
- The child's Y-step reimposes its gangster zeros and box at iteration 1; Z_C is as today.
  The first W = V_C'(Y_C + Z_C/β)V_C then has a positive part close to the parent's
  restricted R. `dnn_bound`'s check that rejects a child start without Y accepts this one;
  the tracker starts with its full ssyevd as today.

The certificate is untouched (Zc still comes from the R-step). With #2, carry the
parent's final β in `DNNWarm`, rescaled by m_C/m.

### 2. Residual-balanced β at refresh points

At each check (every 10 iterations), with L̃ = L/scale:

- r̂ = ‖Y − VRV'‖_F / max(‖Y‖_F, ‖VRV'‖_F);
- ŝ = β‖V'(Y − Y_prev)V‖_F / max(‖V'L̃V‖_F, ‖V'ZV‖_F). The row loop already forms
  ds = Y_new − Y_old; on check iterations write it to a buffer and run one `forward` on it,
  about 3 ms per iteration amortized. The denominators transcribe OSQP's normalization;
- S += log(r̂/ŝ) − log c*, n += 1.

At a scheduled full refresh (every 40 iterations in production), if the iteration is at
most 200, at least 40 iterations have passed since the last change, and fewer than four
changes have been made: f = clip(exp(S/(2n)), 0.5, 2). If f lies outside [1/1.5, 1.5], set
β ← f·β (so ibeta changes), S = n = 0. Z is the unscaled multiplier, so it stays as it is;
only W jumps, and the full ssyevd at the refresh absorbs that, which is why changes wait
for refreshes (OSQP waits for refactorizations for the same reason). The certificate needs
nothing: V'Zc V = β(W − R) ⪯ 0 for the current β. Freezing after 200 iterations gives the
summable-increment condition of Lorenz and Tran-Dinh (2019) and He-Yang-Wang (2000), and
keeps the give-up rule's gap ratio on a stationary iteration.

Target c*: not 1, since the right ratio is scale-dependent and unknown (Wohlberg 2017).
Take c* as the median r̂/ŝ over iterations 100-200 at β = 0.3·m/3 on the 18 all32
fathoms, where that β is the measured optimum. The controller then acts mainly in the
transient, where node 2's ratio is 3x the tail's, so it predicts β about 1.7x higher early
(sqrt 3). The two-phase test in (b) checks this in 10 lines before the controller is built.
Robustness check: it should pick about 0.21 on nug30 children and about 0.3 on wil50; if it
picks the same β on both, it carries no information.

### 3. Scheduled inexactness of the R-step

Inexact ADMM/DR converges when subproblem errors are summable (Eckstein and Bertsekas 1992;
Rontsis et al. 2022's Ritz-residual tolerances) or a fraction of the current fixed-point
residual (Eckstein and Yao 2018); ours is fixed in time. The data fit: one-pass bf16
tracks fp32 to iteration 100 and loses 3.4 at 150 and 7 at 200; refresh 120 and degree 2
show the same late-onset floors [m]. In the node-2 tail the certificate even follows the
refresh phase: 48,607.5-48,608.0 at refresh iterations against 48,603.4-48,605.1 ten
iterations later (npos 1,552 against 1,375-1,537, iterations 1,400-1,460).

Exact change, in `proj_psd_subspace`:

- **Phase E** (from iteration 1, also for warm children): the five W·X products use
  C += Xh·Wh only (one `tdpbf16ps` per tile pair instead of three, and only Wh packed:
  9 MB). CholQR, Rayleigh-Ritz, forming R in fp64 and the refresh are unchanged.
- **Switch:** leave phase E for good when pres_k < 0.25·pres_10 (on node 2 that is about
  iteration 100: pres 0.796 at 10, 0.224 at 100) or at iteration 120, whichever comes
  first. After that, run the 3-pass split as now (error 3.4e-6).

Arithmetic: in bundle A the five products take 81 ms plus a 7 ms pack of 354 ms; one
pass is about a third of the tile work, about 8 ms per product [e], so about 45 ms (13%) is
saved per phase-E iteration: 6% of a 220-iteration cold fathom, 6.5% of a miss given up at
200, less on warm children. Small, but almost free, and the certificate cannot be affected
(R = W − computed negative part; V'Zc V is recomputed in double and Rump-checked).

**Combined:** #1 at its middle (15%) and #3 at 5% give about 1.2x on the dominant node
class; #2 is either a transient gain or a robustness item, and its two-phase test decides.

## (d) References

- Boyd, Parikh, Chu, Peleato, Eckstein 2011, FnT ML 3(1): residual balancing (§3.4.1); over-relaxation 1.5-1.8.
- Stellato, Banjac, Goulart, Bemporad, Boyd 2020, OSQP, Math. Prog. Comp. 12: adaptive ρ (code above), α = 1.6, Ruiz scaling; Table 7 warm start 2.9-3.6x fewer iterations (portfolio).
- Stellato, Naik, Bemporad, Goulart, Boyd 2018, ECC: OSQP in B&B with warm-started children (qualitative).
- O'Donoghue, Chu, Parikh, Boyd 2016, JOTA 169; O'Donoghue 2021, SIAM J. Optim. 31 (SCS 3): adaptive scale (code above), α = 1.5, Anderson reset on rescale.
- Garstka, Cannon, Goulart 2021, COSMO, JOTA 190: chordal decomposition, α = 1.6, Ruiz, adaptive ρ.
- Zheng, Fantuzzi, Papachristodoulou, Goulart, Wynn 2020, CDCS, Math. Prog. 180: chordal ADMM.
- Wen, Goldfarb, Yin 2010, SDPAD, Math. Prog. Comp. 2: SDP ADMM with penalty adjustment.
- Malick, Povh, Rendl, Wiegele 2009, SIAM J. Optim. 20: boundary-point method, σ updated in an outer loop.
- He, Yang, Wang 2000, JOTA 106: self-adaptive penalty, convergence without monotonicity.
- Wohlberg 2017, arXiv 1704.06209: residual balancing is scale-dependent; normalize; target ratio unknown.
- Xu, Figueiredo, Goldstein 2017, AISTATS (PMLR 54): spectral penalty AADMM, incl. DIMACS SDPs (qualitative here).
- Lorenz, Tran-Dinh 2019, Comput. Optim. Appl. 74: non-stationary DR/ADMM, summable step increments.
- He, Liu, Wang, Yuan 2014, SIAM J. Optim. 24 (strictly contractive PRSM); He, Ma, Yuan 2016, SIAM J. Imaging Sci. 9 (symmetric-ADMM domain).
- Graham, Hu, Im, Li, Wolkowicz 2022, INFORMS J. Comput. 34: rPRSM for the QAP DNN.
- Eckstein, Bertsekas 1992, Math. Prog. 55 (relaxed and inexact ADMM as dual DR); Eckstein, Yao 2018, Math. Prog. 170 (relative-error DR/ADMM).
- Rontsis, Goulart, Nakatsukasa 2022, JOTA 192: approximate LOBPCG projections, summable tolerances, up to 20x faster, iterations not noticeably higher.
- Chen, He, Ye, Yuan 2016, Math. Prog. 155: direct multi-block ADMM can diverge.
- Ruiz 2001, RAL-TR-2001-034: row/column equilibration.
