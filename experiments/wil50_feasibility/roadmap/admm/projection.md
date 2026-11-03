# The PSD projection in the DNN ADMM: what the literature outside QAP offers

Research memo, Nov 3 2026; nothing was run. **[m]** = measured (source given),
**[e]** = estimate. Baseline: bundle D of `../../amx/README.md`, 330 ms per
iteration at m = 47, K = 2117. Literature was checked by web search; arxiv.org
could not be fetched from this session, so paper details come from abstracts and
search excerpts and are flagged where only partly confirmed.

## (a) Diagnosis: how far we are from the cheapest defensible projection

**Where the projection's ~245 ms go** (node 10, bundle A, loaded, from
`../../amx/README.md`, with the refresh rescaled to every 40 iterations):

| piece | ms | floor at today's accuracy [e] |
|---|---|---|
| W·Y products, d + 1 = 4 at degree 3, + AMX pack | 88 | 65-75 with per-vector degree (idea 3) |
| one-pass CholQR | 22 | 22 |
| Rayleigh-Ritz: Y'T 13.5, ssyevd(kb ≈ 605) 27.5, X = YG 13 | 54 | 25-35 (idea 2) |
| convert + Chebyshev element-wise | 8 | 8 |
| full ssyevd refresh (734 ms / 40) | 17 | 8-17 |
| form R: dsyrk 35, symmetrize 13, F 7 | 55 | 40 (idea 1) |
| **total** | **~245 (74%)** | **~170-200** |

So the projection is within about **1.25-1.45x** of the cheapest version that
keeps today's accuracy. Over the whole iteration that is 330 → 255-285 ms
(1.15-1.3x). The kernel's unused last third (`../../amx/README.md`) is engineering,
not literature, and is not counted here.

**Today's accuracy is the binding constraint, not a luxury.** Every looser
projection we measured lost iterations: degree 2 left a floor 20 units above
exact at 2,500 iterations; refresh 120 (Weyl bound) lost 81 units at iteration 200;
one-pass bf16 drifted 7 units (`../../../eigensolver/README.md`,
`../../solvers/accel/README.md`, `../../solvers/research_iteration.md`). Degree
3 with tight lambda matches the exact trajectory. Theory says the same thing.
Inexact ADMM/Douglas-Rachford converges when the subproblem errors are
summable (Eckstein-Bertsekas 1992). It also converges under a relative-error
criterion, where the error is bounded by a fixed fraction of the step
(Eckstein-Yao 2018). Our tracker's error after a warm start behaves like
c(d)·||W_k − W_{k−1}||, where c(d) is the filter's damping at the edge of
the spectrum. That is a relative-error rule whose constant is set by the
degree. Degree 2's floor is what a constant that is too large looks like.
Caveats: PRSM with γ1 = γ = 1 is not an averaged operator, so no theorem
covers our exact iteration, and no bound says how small c(d) must be. The
usable lesson: **a cheaper
projection must reach the same edge damping with less arithmetic.** An
"approximate early, exact late" schedule is ruled out by our own data, because
degree 2 was already behind at iteration 50.

**The general-purpose solvers do more work than we do at this size.** SCS and
COSMO compute a full `syevr` of every PSD block on every iteration (as stated
in Rontsis et al. 2022). SDPAD (Wen-Goldfarb-Yin 2010) is a full-eigendecomposition
ADMM, though I could not confirm whether it has a partial-spectrum switch.
SDPNAL+ v1.0 switches to a partial eigendecomposition "whenever it is expected
to be more economical" (Sun-Toh-Yuan-Zhao 2020). That pays only when one side
is small. With a 28%/72% split, all four would run a cold fp64 eigensolve of
order 2117 per iteration. That is about 1.4 s [e], against our ~190 ms warm-started
fp32/bf16 step. The one approximate-projection ADMM with a convergence proof
(Rontsis-Goulart-Nakatsukasa 2022, a modified COSMO) uses warm-started LOBPCG
on the smaller side. It reports up to 20x on SDPLIB with no noticeable increase
in iterations [m]. The paper itself says LOBPCG is efficient only when almost
all eigenvalues share one sign, is slower than an exact eigendecomposition at an
even split, and falls back to a full eigendecomposition otherwise. Our tracker
is already the right family: Chebyshev-filtered subspace iteration across a
sequence of slowly changing matrices (Zhou-Saad-Tiago-Chelikowsky 2006; ChASE,
Winkelmann-Springer-Di Napoli 2019). ChASE's documentation names ≤ 20% of the
spectrum as its sweet spot; we are at 28%.

**Structure gives nothing more.** The 5x10 grid's group has order 4
(`../tree_size.md` §1). A node keeps a symmetry only if every fixed location
is fixed by a group element. For the row reflection that means every fixed
location is in row 2. That is about (10/50)³ ≈ 1% of depth-3 nodes, and even
then the group has order 2: two blocks, halving the products (de
Klerk-Sotirov 2010). Chordal decomposition needs a sparse aggregate pattern.
Ours is dense, because Z is dense and Y ≥ 0 couples every entry. No choice of
basis for V block-diagonalizes W; the R-step is invariant under the basis
anyway (`dnn.C` header §1).

## (b) Ranked ideas

| # | idea | effect per iteration | days | risk | evidence | cheapest first test → kill |
|---|---|---|---|---|---|---|
| 1 | **fp32 rank term on non-certificate iterations only**: R = W + fl32(FF') on 9 of 10 iterations, fp64 dsyrk when `rstep_cert` | −12 to −13 ms (4%) | 0.5 | low | [m] `DNN_R32=1` everywhere: 409 → 395 ms, same fathom; it was closed only because of the certificate margin (5e-6 → 0.08), which this removes | node 10, 300 iterations, and depth-4 node 3. **Kill:** certificate margin above 1e-5, a fathom moves, or the bound at 300 moves by more than 0.1 |
| 2 | **Complementary-subspace Rayleigh-Ritz** (CS2, Banerjee et al. 2018): P_−(H) = H − P_+(H), computing only the kb − nneg ≈ 60-110 positive Ritz pairs of the 605² matrix H; F from a pivoted Cholesky; the next block is the CholQR basis, with no rotation | −15 to −30 ms (5-9%) [e] | 2-3 | low-med | Same R as today in exact arithmetic; CS2 cuts subspace diagonalization to the few top states (DFT, JCTC 2018; no headline factor found) | Microbenchmark on this VM, pinned and repeated: ssyevd(605) + the K × 605 × 605 rotation + scaling vs ssyevr(605, positive part) + spstrf + trapezoidal F. **Kill:** saving below 10 ms. Then node 10 / all32: **kill** on any fathom change |
| 3 | **Per-vector Chebyshev degree** (ChASE's degree optimization): strongly negative Ritz vectors get degree 1-2, the edge keeps degree 3 | −10 to −25 ms (3-8%) [e] | 1-2 | medium | ChASE reports its degree optimization cuts FLOPs (figure not checked); our error weighting argument (§c3) | Free: log t_i = (θ_i − c)/h per iteration on node 10 and compute the planned column count. **Kill:** fewer than 25% of the block's columns drop a degree. Then all32: **kill** on any fathom later, or misses' lb@200 lower by more than 1 |
| 4 | Free residual monitor ‖E‖²_F = ‖WY‖²_F − Σθ_i² (Y orthonormal), with refresh on demand instead of every 40 iterations; gap-free error bound of Goulart-Nakatsukasa-Rontsis 2020 | −5 to −12 ms (refresh 17) [e] | 1-2 | medium | The monitor is O(K·kb) and exact algebra; it does not see negative eigenvalues outside the block, which is the refresh's job | Log the monitor next to the refresh's correction on node 2 to 2,500 iterations. **Kill:** refreshes still needed more often than every 60 iterations |
| 5 | Accuracy schedule from inexact-ADMM theory (loose early, tight late) | iterations lost [m] | – | – | Degree 2 with tight λ was behind at iteration 50, 100, 150 and 200 | Dead |
| 6 | LOBPCG / Davidson, warm-started (Rontsis et al., COSMO prototype) | slower [e] | – | – | Rayleigh-Ritz on [X, R, P] = 1,800 columns ≈ 0.6 × a full 2117 eig, about 450 ms; LOBPCG's own authors exclude near-even splits | Dead by arithmetic |
| 7 | FEAST / contour filters | ~190 GF vs our ~22 GF [e] | – | – | Each quadrature node is a dense complex LU (~25 GF) plus a K×600 solve; 4-8 nodes | Dead: rational filters pay only with sparse solves |
| 8 | Randomized range finder / Nyström sketch (Halko-Martinsson-Tropp 2011) | none | – | – | This is subspace iteration from a random start: our method without the warm start | Dead |
| 9 | Burer-Monteiro, SDPLR, factored or low-rank ADMM; spectral bundle | none | – | – | rank R* ≈ 1,500 of 2,117; Y ≥ 0 cannot be factored; a fast bundle needs size ≥ rank of the optimal primal (Ding-Grimmer 2023). Our negative form already keeps the dual's PSD part factored (F, K × nneg) | Settled in `../../solvers/research_solvers.md` |
| 10 | Composite polynomial / Newton-Schulz / polar sign iterations | 420-590 GF [e] | – | – | | Settled in `../../solvers/research_iteration.md` item 4 |

Ideas 1-3 together: about −40 to −65 ms, 330 → 265-290 ms (1.15-1.25x) [e].
They touch different lines of `proj_psd_subspace` and stack.

## (c) The top three in detail

### 1. fp32 R except on certificate iterations

The certificate uses only the R-step of its own iteration:
`Zc = Z + β(Y − VRV')`, so V'ZcV = β(W − R). It does not depend on any
earlier R. So R needs W − R = −fl64(FF') only when `rstep_cert` is true.
Everywhere else, R is just the ADMM step. A rounding of 1e-7 relative is
orders below the tracker's own error: TF32-level operands track fp32 exactly
(`research_iteration.md` §2).

**Change.** Add a `bool cert_step` argument to `proj_psd_subspace`, passed as
`rstep_cert` from the main loop. Replace `if (r32 && k > 0)` with
`if ((r32 == 1 || (r32 == 2 && !cert_step)) && k > 0)` and make
`DNN_R32=2` the production value. Nothing else changes: the fp32 path already
fuses the symmetrization into the write.

The certificate and the double re-projection path are untouched; certificate
iterations run today's fp64 path, so their margin must stay at about 5e-6.

### 2. Complementary-subspace Rayleigh-Ritz

After CholQR, Y (K × kb) is orthonormal, T = WY, and H = Y'T. Today we
eigendecompose all of H (27.5 ms), rotate X = YG (13 ms), and scale F (7 ms).
But R only needs

R = W − Y P_−(H) Y' = W + Y S Y',  where S = P_+(H) − H ⪰ 0,

and P_+(H) involves only the npos_H = kb − nneg positive Ritz pairs. These
are the guard band, about 60-110 of 605. This is the CS2 identity (Banerjee
et al. 2018, Kohn-Sham density matrices). There it removes the full subspace
diagonalization and keeps only the few top states.

**Change**, in `proj_psd_subspace` after the `G = Y'T` symmetrization:

1. Copy G to H2. Call `ssyevr("V","V","L", kb, H2, vl = 0, vu = FLT_MAX)`.
   This returns the m = npos_H positive pairs U_+, θ_+. Then nneg = kb − m and
   `theta_top = θ_+[m−1]`.
2. S = U_+ diag(θ_+) U_+' − G: an ssyrk of rank m on sqrt(θ_+)-scaled U_+,
   minus G. Cost kb²·m ≈ 40 MF.
3. Call `spstrf` on S with tolerance 1e-6·max|θ|. This gives rank r ≈ nneg
   and pivots P. Set F = Y(:,P)·L, K × kb × r and trapezoidal (strmm on the
   leading r × r block plus sgemm for the rest). F replaces both the X = YG
   update and the scaling loop.
4. R = W + FF' as today: dsyrk on certificate iterations, ssyrk otherwise
   (idea 1). W − R = −FF' is negative semidefinite by construction whatever
   the accuracy of U_+ or L. **The certificate is unaffected.**
5. The next iteration filters Y (pointer swap into `w.Xs`), not the Ritz
   vectors. Subspace iteration is basis-independent, and CholQR (orthogonal
   iteration) keeps Y close to ascending Ritz order, so the next CholQR
   stays well conditioned.
6. Block growth works unchanged: append random vectors. When the block must
   shrink (kb > nneg + 1.5·guard), or every 10th iteration as a safeguard,
   run today's full Rayleigh-Ritz path. That path re-sorts the basis and trims it.

Estimate: 47.5 ms today → 20-34 ms [e]. ssytrd is still in ssyevr, and at
n = 605 it is half Level-2, which is why the range is wide. Once this works,
the second stage of CS2 removes ssytrd: a small Chebyshev filter on H for its
top m eigenpairs. It starts from the last m unit vectors (Y is near Ritz
order) and costs 4-6 products of 605 × 605 × 110, about 5 ms. In exact
arithmetic the iterate is the same as today's, so the 32-node run should
reproduce bundle D's decisions exactly.

### 3. Per-vector filter degree

The filter damps [theta_top, lam_hi] by 1/|C_d(t)|, with t_i = (θ_i − c)/h.
Column i of the block contributes about |θ_i|·sin φ_i to the error in R.
After a warm start, its leakage into the positive eigenspace is about
‖ΔW‖/(|θ_i| + λ⁺). The error contribution after filtering is therefore about
‖ΔW‖/|C_d(t_i)|, which falls steeply as θ_i becomes more negative. Strongly
negative columns reach the edge column's damping at a lower degree. As an
example, take theta_top = 0.1 and lam_hi = 1.3, so c = 0.7 and h = 0.6.
Degree 3 at the edge (θ ≈ 0) damps by 2.9. Degree 2 gives at least 2.9 for
θ ≤ −0.14, and degree 1 for θ ≤ −1.04. ChASE (Winkelmann et al. 2019)
optimizes per-vector degrees in the same way for sequences of eigenproblems.

**Change**, in `proj_psd_subspace` before the filter loop. The block X is in
ascending Ritz order with values `theta[]` from the previous Rayleigh-Ritz.
Under idea 2, use diag(H) of the previous iteration instead.

- Set D = `degree` and t_e = t of the largest negative θ. For each column,
  d_i = min{d ≥ 1 : |C_d(t_i)| ≥ |C_D(t_e)|}. Guard and positive columns get
  D. Because d_i does not decrease with i, step d is active on a trailing
  range [j_d, kb).
- The first product runs on all kb columns. For d = 2..D, call
  `WY(Y + j_d*K, T + j_d*K)` with kb − j_d columns (`amx_gemm` takes any
  column count; the sgemm call takes `kb - j_d`). Run the recurrence on that
  range only. When a range retires, copy its columns into the buffer that will
  hold the final Y. The Y/Yp swaps leave them in alternating buffers.
- Before CholQR, scale every column to unit norm (O(K·kb)). Columns filtered
  to different degrees differ in norm, and one-pass CholQR should not see that.

Products then cost kb + Σ_{d≥2} (kb − j_d) columns plus the Rayleigh-Ritz
product, instead of (D + 1)·kb. If a third of the negative block retires at
degree 1 and another third at degree 2, the products drop from 4 kb to 3 kb:
−20 ms [e]. The histogram test in the table settles this before any timing
run. The risk is iteration count, not validity: the certificate never sees
the filter.

## (d) References

- Eckstein & Bertsekas, Math. Program. 55:293-318, 1992. DR splitting and ADMM
  as proximal point methods; generalized ADMM with errors in the subproblems
  (summable-error theorem in §3, not re-read here).
- Eckstein & Yao, Math. Program. 170:417-444, 2018. Relative-error approximate
  DR and ADMM.
- Goulart, Nakatsukasa & Rontsis, Linear Algebra Appl. 594:177-192, 2020.
  Gap-independent bound on the error of an approximate PSD projection from
  Ritz residuals; the exact constant was not confirmed (excerpt suggests √2‖R‖).
- Rontsis, Goulart & Nakatsukasa, JOTA 192:292-320, 2022 (arXiv 1912.02767).
  Approximate ADMM in COSMO with warm-started LOBPCG; up to 20x on SDPLIB, no
  noticeable increase in iterations [m]; worst at an even eigenvalue split.
- Garstka, Cannon & Goulart, JOTA 190:779-810, 2021 (COSMO): full `syevr` per
  iteration. O'Donoghue, Chu, Parikh & Boyd, JOTA 169, 2016 (SCS): full
  eigendecomposition per iteration.
- Wen, Goldfarb & Yin, Math. Program. Comput. 2:203-230, 2010 (SDPAD).
- Sun, Toh, Yuan & Zhao, Optim. Methods Softw. 35, 2020 (SDPNAL+ v1.0):
  partial eigendecomposition when economical.
- Souto, Garcia & Veiga, Optimization 71(1), 2022 (ProxSDP): low-rank
  approximate projection in a PDHG variant, with convergence guarantees.
- Zhou, Saad, Tiago & Chelikowsky, J. Comput. Phys. 219:172-184, 2006.
  Chebyshev-filtered subspace iteration across SCF steps; a filter step 5-10x
  cheaper than an eigensolver at the same step count [m].
- Winkelmann, Springer & Di Napoli, ACM TOMS 45(2), 2019 (ChASE, arXiv
  1805.10121). Per-vector degree optimization; efficient up to ~20% of the
  spectrum.
- Banerjee, Lin, Suryanarayana, Yang & Pask, J. Chem. Theory Comput.
  14:2930-2946, 2018 (CS2). Complementary subspace: diagonalize only the top
  states, with an inner Chebyshev filter.
- Polizzi, Phys. Rev. B 79:115112, 2009 (FEAST). Knyazev, SIAM J. Sci. Comput.
  23, 2001 (LOBPCG). Halko, Martinsson & Tropp, SIAM Rev. 53, 2011.
- Ding & Grimmer, SIAM J. Optim. 33:1305-1332, 2023. Spectral bundle methods
  are linear only when the bundle captures the optimal rank.
- de Klerk & Sotirov, Math. Program. 122:225-246, 2010. Symmetry reduction of
  QAP SDP relaxations.
- In-repo: `latest/dnn.C`, and the notes linked above.
