# Certificate study: the arguments (companion to REPORT.md)

Notation as in `docs/dnn_bound/dnn_bound.tex` section 2: lifted y = (1, x),
face basis V, Zc = Z + beta (Y - VRV') the R-step multiplier (PRSM half step),
M = L + scale Zc (plus L_add for cuts), and the certificate

    lb = s + GL(M) - (m+1) scale lambda - margin,   lambda >= lambda_max(V'ZcV),

GL(M) = M00 + LAP_m(c), c_ij = M_0,ij + M_ij,0 + M_ij,ij + LAP_{m-1}(block row ij).

## 1. Where the certificate can lose, at a fixed iterate

1. The semidefinite term (m+1) scale lambda and the rounding margins. In every
   wil50 log we have they total 3.6e-6 to 4.95e-6 units, with no failed
   Cholesky test and no refinement (`amx/logs/val_*`, `warm_y/runs/child_9.log`).
   V'ZcV = beta (W - P(W)) is negative semidefinite by construction, so lambda
   is a rounding-level number. Nothing to gain from sharper eigenvalue bounds,
   Rump-type or otherwise.
2. The inner LAPs: the row-column reduction is used until the value is within
   1e-3 |cutoff| (48.8) of the cutoff, then exact Hungarian solves (AUTO). The
   fathoming checks are therefore already exact; reduce-mode values only feed
   the give-up rule and the Polyak targets of the cut phase.
3. The polyhedral relaxation itself: GL is the minimum of <M, Y> over P_GL
   (column blocks x_ij times doubly stochastic, x doubly stochastic), which does
   not see that Y is symmetric. This is the only place with real slack at a
   non-converged iterate (section 3).
4. The multiplier Zc itself (dual not converged): only more iterations, or a
   different multiplier, change this (sections 4, 5).

## 2. The free dual variables are already optimal or reduce to a skew matrix

For every lifted permutation y, <N, yy'> = 0 for three families, so adding any
of them to M keeps the certificate valid (the semidefinite part is untouched):

- (a) the shift t I: <I, yy'> = m+1. Adding t I to Zc raises GL by exactly
  (m+1) t and lambda by t, so the certificate is unchanged. Neutral.
- (b) the Y00 multiplier t E00 (the variable Kim-Kojima-Toh bisect on in the
  Lagrangian-DNN bound): GL rises by t, but V'(Zc + tE00)V = V'ZcV + t vv'
  with v = V'e0, and V'ZcV is exactly zero on range(P(W)) (the complement of
  its negative eigenspace), which contains a part of v of squared norm about
  1/2 (Y00 = 1 sits in R's range). So lambda >= t ||P v||^2 and the penalty
  (m+1) t ||Pv||^2 exceeds t for m >= 2. No gain at any iterate of this method.
- (c) the assignment-constraint multipliers N = A'G + G'A (A: the row and
  column sums of x equal y0). With exact LAPs GL is invariant under A'G (the
  constraints applied to the columns GL reads: GL's own LAP duals absorb
  them), and G'A = A'G + (G'A - A'G) differs from it by a skew-symmetric
  matrix. So the whole freedom is max over skew K of GL_exact(M + K), the
  Hahn-Grant / Adams-Johnson RLT1 dual: K is valid because <K, yy'> = 0 for
  every symmetric yy'.

So "a better choice of the free dual variables by a small exact problem" is
exactly the skew polish (bits 2, 32 in `latest/dnn.C`); there is no separate
LP over shift or assignment multipliers to solve. Its ceiling is the RLT1 LP
min over P_GL ∩ {Y = Y'} of <M, Y>, which is at most the DNN value (the DNN
feasible set lies in it), so it gains only where Zc is not dual optimal.

### The polish (hg_polish) and why every value is a valid bound

R starts as M. A pass: (1) for each block row (ij), an LAP (exact or
row-column reduction) with duals u, v leaves R_(ij),(kl) - u_k - v_l >= 0 and
moves sum u + sum v to the linear term of (ij): an identity for permutations,
because a permutation with x_ij = 1 uses each k != i once and each l != j once
in that row; (2) the outer LAP on c with duals U, V moves sum U + sum V to the
constant B and leaves c - U_i - V_j >= 0: an identity because x is a
permutation; (3) the leftover c'_ij is spread as c'_ij / (m-1) over block row
(ij): an identity because that row has exactly m-1 ones when x_ij = 1;
(4) every non-gangster pair R_ab, R_ba is replaced by its average: an
identity because yy' is symmetric. After any number of passes,
<M, yy'> = B + <R, yy'> >= B for every lifted permutation, since R >= 0 where
yy' can be 1. So B - (m+1) scale lambda - margins is a certificate. Rounding:
the first pass is exactly the exact GL of M (certify's margin formula applies
with its magnitudes); later passes add an explicit allowance (2u times the
magnitudes of the m(m-1) block entries and m linear terms a permutation reads,
times the number of roundings, doubled). Checked: random n = 8 enumeration,
five 12-size QAPLIB instances at every check (tai12b: the polish lifts the
converged certificate by 615 and stays 4 below the optimum).

The skew part of the transfers can be accumulated in an explicit K (exactly
skew: K_ba = -K_ab). With exact LAPs GL_exact(M + K) >= B (the reductions and
the spreading are face terms that GL's own LAPs absorb); measured +0.17 to
+7.1 on nug20. Carrying K to the next check (bit 64) is valid but makes the
bound worse than the plain exact GL (-24 at iteration 100 on nug20): a K tuned
to one iterate's noise is noise at the next.

## 3. Other multipliers from the same iterate

- Line search Zc + t (Z_{k+1} - Zc) (bit 4). Z_{k+1} is exactly box-dual
  feasible (L + Z_{k+1} is zero on the entries the Y-step leaves strictly
  inside [0, 1]) but not semidefinite-feasible; Zc the reverse. The certificate
  g(t) is concave in t (GL concave, lambda_max convex), so a line search
  needs no safeguard; lambda by Lanczos (40 steps) and Rump's test at a shift
  above the Ritz value. Measured on nug12, nug20 (and E1): g decreases from
  t = 0 at every check. Dead.
- Averages of the last k multipliers (bit 8): g(mean Zc_i) >= mean g(Zc_i) by
  concavity, and V'(mean Zc)V <= 0, but while the certificate rises the
  average trails the last iterate (nug20: -8 to -130). Tiny gains (+0.0006 on
  nug12) only at convergence, where nothing is decided. Dead.

## 4. Child bounds from the parent's multiplier (question 3)

For a child fixing facility a to location b, every child permutation is a
parent permutation, so for the parent's Zc and any P containing the parent's
lifted permutations,

    child optimum >= min over P ∩ {x_ab = 1} of <M, Y> - (m+1) scale lambda.

Today's child_lb is GL + r_ab (outer reduced cost); fixing the outer row
exactly gives GL(M | outer x_ab = 1) >= it, and fixing a, b in the inner
LAPs too (scripts/child_gl.py: c^ab_kl = lin_kl + M_kl,ab + M_ab,kl +
LAP_{m-2}) gives the exact per-child GL, never lower. All are valid by the
same proposition as the parent's certificate.

**The ceiling is complementary slackness.** At a dual-optimal Z*, GL(L + Z*)
equals the DNN value, attained by the DNN optimum Y* in P_GL; the minimizing
face of P_GL contains Y*, so it contains a vertex with x_ab = 1 for every cell
with x*_ab > 0, and that child's restricted minimum equals the parent's. The
same holds for the RLT1 (skew-polished) polytope. So child bounds from a
parent's near-converged multiplier exceed the parent's value only on cells
the DNN already puts (almost) no mass on; those children are the cheapest to
fathom anyway. Away from convergence the child bounds are less flat, but the
parent's own value is lower by more.
